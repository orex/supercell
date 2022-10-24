#!/usr/bin/python3
import itertools
import re
import numpy as np
import scipy.stats as stats
from collections import OrderedDict
import argparse
import pathlib
import tempfile
import os
import subprocess
import copy
from pymatgen.core import Structure
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

general_format='-i {cif} -s {cell} -v 1 -q -m -n l{smp} -n h{smp} -n r{smp}'

run_cases = [
    {"format": general_format + " -g",
     "params": OrderedDict([
         ("cif", ["RB-PST-1-DEHY_1_new.cif"]),
         ("cell", ["1x1x1"]),
         ("smp", ["1", "30", "100", "1000"])]),
     "checks": {"es", "syms", "rand"}
     },
    {"format": general_format + " -g",
     "params": OrderedDict([
         ("cif", ["SrSiAlO_1x2x1.cif"]),
         ("cell", ["1x1x1"]),
         ("smp", ["100"])]),
     "checks": {"sort", "es", "syms", "rand"}
     },
    {"format": general_format + " --random-seed {seed}",
     "params": OrderedDict([
         ("cif", ["PbSnTe2.cif"]),
         ("cell", ["2x2x2"]),
         ("seed", ["1546733744", "538425541", "92928180", "2074289474"]),
         ("smp", ["300", "3000"])]),
     "checks": {"syms", "rand"}
     },
    #  Very long
    # {"format": general_format,
    #  "params": OrderedDict([
    #      ("cif", ["CaAl6Te10.cif"]),
    #      ("cell", ["1x1x1"]),
    #      ("smp", ["2000"])]),
    #  "checks": {"es", "syms", "rand"}
    #  }
]

def run_supercell(cmd, params, outdir):
    fullcmd="{cmd} {params} -o {out}".format(cmd=cmd, params=params, out=outdir / "electro")
    print(fullcmd)
    result = subprocess.run(fullcmd, capture_output=True, text=True, shell=True)

    if result.returncode != 0:
        print(result.stderr)
        return None, None

    tot_syms = re.search("([0-9]+) +symmetry operation found for supercell.*", result.stdout)
    tot_comb = re.search("Combinations after merge: ([0-9]+)", result.stdout)
    return outdir, int(tot_syms.group(1)), int(tot_comb.group(1))

def check_case(checks, cmd, fmt, params, data_folder, tmpdir):
    fosuffix = "test_electro"+fmt.format(**params).replace(" ", "_").replace(".", "_").replace("-", "")
    workdir = tmpdir / fosuffix
    os.mkdir(workdir.resolve())
    p=copy.deepcopy(params)
    p["cif"]=(data_folder / p["cif"]).resolve()
    outdir, tot_syms, tot_comb=run_supercell(cmd, fmt.format(**p), workdir)

    cm = {}
    result=True
    for fn in os.listdir(outdir):
        m = re.match("electro_i([lhr])([0-9]+)_w([0-9]+).*\.cif", fn)
        if m:
            # print(fn)
            with open(outdir / fn, "r") as f:
                data = f.read()
                mx = re.match("# E_e = ([-0-9\.]+) eV",data.split("\n")[1])
                se = float(mx.group(1))
                structure = Structure.from_str(data, fmt="cif")
                f.close()

            cm.setdefault(m.group(1), []).append(int(m.group(2)))
            if "es" in checks:
                ew = EwaldSummation(structure, acc_factor=5.0)
                t = ew.total_energy
                # The difference is due to coordinate rounding in cif file.
                if abs(se - t) > 1E-2:
                    print("FAILED energy match('{}'): {} vs {}".format(fn, t, se))
                    result=False

            if "syms" in checks:
                str_syms = len(SpacegroupAnalyzer(structure).get_space_group_operations())
                if tot_syms / str_syms != int(m.group(3)):
                    print("Symmetry operation problems: {}".format(fn))
                    result=False

    if "sort" in checks:
        tot_e = []
        with open(outdir / "electro_coulomb_energy.txt", "r") as f:
            for l in f.readlines():
                m = re.match(".*/electro_i([0-9]+).*\.cif\s+([0-9\.\-]+)\s+eV", l)
                if m:
                    tot_e.append((int(m.group(1)), float(m.group(2))))
            f.close()
        tot_e.sort(key=lambda x: x[1])
        cm["l"].sort()
        cm["h"].sort()
        if cm["l"] != sorted([x[0] for x in tot_e[0:len(cm["l"])]]) or cm["h"] != sorted(
                [x[0] for x in tot_e[-len(cm["h"]):]]):
            print("Bad energies sorting")
            result = False

    if "rand" in checks:
        if len(cm["r"]) != int(p["smp"]):
            print("Failed number of random structures")
        if stats.kstest(np.array(cm["r"]) / (tot_comb - 1), stats.uniform.cdf).pvalue < 0.03:
            print("Bad randomness of structures")

    return result

def main(args):
    tmpdir = pathlib.Path(tempfile.mkdtemp())

    for c in run_cases:
        for t in itertools.product(*[l for l in reversed(c["params"].values())]):
            fv = OrderedDict([(k, v) for k, v in zip(reversed(c["params"].keys()), t)])
            check_case(c["checks"], args.supercell_cmd, c["format"], fv, args.data_folder, tmpdir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create tests.')
    parser.add_argument('--supercell-cmd', type=str, default='supercell', help='Path to checked supercell program')
    parser.add_argument('--data-folder', type=pathlib.Path, default=pathlib.Path("./"), help='Path with data files.')
    parser.add_argument('--tmp-folder', type=pathlib.Path, default=pathlib.Path(tempfile.gettempdir()), help='Path to store temporary files')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    print(args)
    if args:
        main(args)
