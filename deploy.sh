#!/bin/bash

function deploy-doc {
  cp ${c_path}/build/doc/man/supercell_man.{pdf,html} ${REPO_PATH}/${d_prefix}/doc/.
  cp ${c_path}/doc/man/supercell_man.css ${REPO_PATH}/${d_prefix}/doc/.
  cp ${c_path}/build/doc/tutorial/supercell_tutorial.pdf ${REPO_PATH}/${d_prefix}/doc/.
}

function deploy-exe {
  tmp_folder=`mktemp -d -t XXXXXX`
  cd ${tmp_folder}
  cp ${c_path}/build/src/sc_cli/supercell* .

  wget -nv https://github.com/orex/supercell/raw/deploy/README -O README

  if [[ "$1" == "windows" ]]; then
    zip -9 supercell-$1.zip *
  else
    tar czvf supercell-$1.tar.gz *
  fi

  cp supercell-$1.* ${REPO_PATH}/${d_prefix}/exe/.
}


set -e # Exit with nonzero exit code if anything fails

# Pull requests and commits to other branches shouldn't try to deploy, just build to verify
if [ "$DEPLOY_BUILD" == "FALSE" -o "$DEPLOY_BUILD" == "" ]; then
    echo "Skipping deploy."
    exit 0
fi

# Save some useful information
DEPLOY_REPO="git@github.com:orex/orex.github.io.git"
DEPLOY_BRANCH="master"
DEPLOY_DIR="supercell"
COMMIT_SHA=`git rev-parse --verify HEAD`

ssh-keyscan -t rsa github.com > ~/.ssh/known_hosts
echo "${DEPLOY_KEY}" > ~/.ssh/id_rsa 2>/dev/null

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deply)
git clone --recurse-submodules "${DEPLOY_REPO}" --branch "${DEPLOY_BRANCH}" out && cd out
REPO_PATH=${PWD}
d_prefix="${DEPLOY_DIR}"
if [ "$TRAVIS_BRANCH" != "master" ]; then
  d_prefix="${DEPLOY_DIR}/deploy_artifacts/$TRAVIS_BRANCH"
fi

# Now let's go have some fun with the cloned repo
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

#Main code
mkdir -p ${d_prefix}/doc
mkdir -p ${d_prefix}/exe

if [[ "$DEPLOY_BUILD" == "doc" ]]; then
  deploy-doc
else
  deploy-exe "$DEPLOY_BUILD"
fi

cd ${REPO_PATH}

ls -lah

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add --no-ignore-removal .

git status

git commit -m "Deploy to GitHub Pages $TRAVIS_OS_NAME: ${COMMIT_SHA}"

# Now that we're all set up, we can push.
git push $DEPLOY_REPO $DEPLOY_BRANCH
