#!/usr/bin/python3
import itertools
import re
from collections import OrderedDict
import argparse
import pathlib
import tempfile
import os
import random
import subprocess
import copy
import filecmp

general_format='-i {cif} -s {cell} -v 1 --random-seed {seed} -q -m -n f{smp} -n a{smp} -n l{smp} -n h{smp} -n w{wstore}'

match_lines_re = [
    "Random SEED:.*([0-9]+).*",
    "Total charge cell: +([-0-9\.]+).*",
    "Chemical formula of the supercell: +(.*)",
    "Total charge of supercell: +([-0-9\.]+).*",
    "Minimal distance between atoms of two distinct groups: +([0-9\.]+) +A.",
    "The total number of combinations is +([0-9]+).*",
    "([0-9]+) +symmetry operation found for supercell.*",
    "Combinations after merge: ([0-9]+)"
]

def full_cpu_check():
    ncpusstr = subprocess.run('lscpu -p=CPU', capture_output=True, text=True, shell=True).stdout.split('\n')
    ncpus = [int(d) for d in ncpusstr if d.isnumeric()]
    if len(ncpus) == 1:
        return ncpus
    if ncpus == 2:
        return [[ncpus[random.randint(0, 1)]], ncpus]
    ri = copy.deepcopy(ncpus)
    random.shuffle(ri)
    ri=ri[0:random.randint(2, len(ncpus) - 1)]
    ri.sort()
    return [[ncpus[random.randint(0, len(ncpus) - 1)]], ri, ncpus]

run_cases = [
    ({"format" : general_format, "params": OrderedDict([
        ("cif"    , ["Ca2Al2SiO7.cif", "PbSnTe2.cif"]),
        ("cell"   , ["1x1x1", "1x1x2", "1x2x1", "1x2x2", "2x1x1", "2x1x2", "2x2x1", "2x2x2"]),
        ("smp"    , ["1", "10", "256", "1000", "10000"]),
        ("seed",    ["2087334979", "1893125674", "1342325495", "548259826", "1027322602"]),
        ("wstore" , ["0"])
        ]
    ), "cpucheck": full_cpu_check}),
    ({"format" : general_format, "params": OrderedDict([
        ("cif"    , ["RB-PST-1-DEHY_1_new.cif"]),
        ("cell"   , ["1x1x1"]),
        ("smp"    , ["1", "10", "256", "1000", "10000"]),
        ("seed",    ["2087334979", "1893125674"]),
        ("wstore" , ["0"])
    ]
    ), "cpucheck": lambda: []}),
]

def run_supercell(cpumask, cmd, params, id, outdir):
    sfx="" if len(cpumask) == 0 else "taskset -c {} ".format(cpus_to_str(cpumask, ','))

    outfname=outdir / "{id}-mask{cpumask}.tar".format(id=id, cpumask=cpus_to_str(cpumask, '-'))
    fullcmd="{sfx} {cmd} {params} -a {tar}".format(sfx=sfx, cmd=cmd, params=params, tar=outfname)
    result = subprocess.run(fullcmd, capture_output=True, text=True, shell=True)

    if result.returncode != 0:
       print(result.stderr)
       return None, None

    return outfname, result.stdout

def cpus_to_str(cpus, join):
    return "" if cpus is None or len(cpus) == 0 else join.join([str(i) for i in cpus])

def check_case(refcmd, testcmd, cpumasks, fmt, params, data_folder, tmpdir):
    fosuffix = "test_"+fmt.format(**params).replace(" ", "_").replace(".", "_").replace("-", "")
    workdir = tmpdir / fosuffix
    os.mkdir(workdir.resolve())
    p=copy.deepcopy(params)
    p["cif"]=(data_folder / p["cif"]).resolve()
    refinput=run_supercell([], refcmd, fmt.format(**p), "ref", workdir)
    all_input=OrderedDict([("ref", refinput)])
    for m in cpumasks:
        all_input["test-mask{}".format(cpus_to_str(m, '-'))] = run_supercell(m, testcmd, fmt.format(**p), "test", workdir)

    rest={}
    for l in all_input.values():
        axs = l[1].split('\n')
        for i, r in enumerate(match_lines_re):
           rest.setdefault(i, []).append([re.match(r, m).group(1) for m in axs if re.match(r, m)])

    result=True
    for k, v in rest.items():
        if v.count(v[0]) != len(v):
            result=False
            print("Problems matching '{}'".format(match_lines_re[k]))

    for k, v in all_input.items():
        if k != "ref" and not filecmp.cmp(workdir / all_input["ref"][0], workdir / v[0], shallow=False):
            print("Problems output '{}'".format(v[0]))
            result=False

    return result


def main(args):
    index=0
    if args.action == 'test':
        tmpdir = pathlib.Path(tempfile.mkdtemp())

    for c in run_cases:
        for t in itertools.product(*[l for l in c["params"].values()]):
            if len(args.testindexes) == 0 or index in args.test_numbers:
                fv = OrderedDict([(k, v) for k, v in zip(c["params"].keys(), t)])
                passed=True
                if args.action == 'test':
                    passed = check_case(args.supercell_ref, args.supercell_test, c["cpucheck"](), c["format"], fv,
                               args.data_folder, tmpdir)
                if args.action == 'list' or args.verbose or not passed:
                    print(("{index}: supercell {params} {fail}").format(index=index, params=c["format"].format(**fv),
                          fail=("" if passed else ": FAILED!")))

        index += 1



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create tests.')
    parser.add_argument('testindexes', metavar=int, type=int, nargs='*',
                        help='Test numbers to run. Default: all')
    parser.add_argument('-a', '--action', choices=['list', 'test'], default='list', help='Action to perform')
    parser.add_argument('--supercell-ref', type=str, default='supercell', help='Path to reference supercell program')
    parser.add_argument('--supercell-test', type=str, required=True, help='Path to checked supercell program')
    parser.add_argument('--data-folder', type=pathlib.Path, default=pathlib.Path("./"), help='Path with data files.')
    parser.add_argument('--tmp-folder', type=pathlib.Path, default=pathlib.Path(tempfile.gettempdir()), help='Path to store temporary files')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    print(args)
    if args:
       main(args)
