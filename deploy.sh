#!/bin/bash

function deploy-doc {
  cp ${c_path}/build/doc/man/supercell_man.pdf ${deploy_dir}/doc/.

  cp ${c_path}/build/doc/man/supercell_man.html ${deploy_dir}/doc/.
  cp ${c_path}/doc/man/supercell_man.css ${deploy_dir}/doc/.

  cp ${c_path}/build/doc/tutorial/supercell_tutorial.pdf ${deploy_dir}/doc/.

}

function deploy-exe {
  tmp_folder=`mktemp -d -t XXXXXX`
  cd ${tmp_folder}
  cp ${c_path}/build/src/sc_cli/supercell .

  for i in {atomtyp.txt,bondtyp.txt,element.txt,phmodel.txt,space-groups.txt,types.txt}
  do
    wget -nv https://github.com/openbabel/openbabel/raw/master/data/$i -O $i
  done

  wget -nv https://github.com/orex/supercell/raw/deploy/README -O README

  tar czvf supercell-$1.tar.gz *

  cp supercell-$1.tar.gz ${deploy_dir}/exe/.
}


set -e # Exit with nonzero exit code if anything fails

SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"

# Pull requests and commits to other branches shouldn't try to deploy, just build to verify
if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" -o "$DEPLOY_BUILD" != "TRUE" ]; then
    echo "Skipping deploy."
    exit 0
fi

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deply)
git clone $REPO out
cd out
deploy_dir=${PWD}
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH

# Now let's go have some fun with the cloned repo
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"


#Main code

mkdir -p doc
mkdir -p exe

wget -nv https://github.com/orex/supercell/raw/deploy/index.html -O index.html

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  deploy-doc
  deploy-exe linux
elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
  deploy-exe osx
fi

cd ${deploy_dir}

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add --no-ignore-removal .
git commit -m "Deploy to GitHub Pages $TRAVIS_OS_NAME: ${SHA}"

# Get the deploy key by using Travis's stored variables to decrypt deploy_key.enc
ENCRYPTION_LABEL="afb0f7bfe5ee"
ENCRYPTED_KEY_VAR="encrypted_${ENCRYPTION_LABEL}_key"
ENCRYPTED_IV_VAR="encrypted_${ENCRYPTION_LABEL}_iv"
ENCRYPTED_KEY=${!ENCRYPTED_KEY_VAR}
ENCRYPTED_IV=${!ENCRYPTED_IV_VAR}

wget -nv https://github.com/orex/supercell/raw/deploy/id_rsa_gh-pages.enc
openssl aes-256-cbc -K $ENCRYPTED_KEY -iv $ENCRYPTED_IV -in id_rsa_gh-pages.enc -out id_rsa_gh-pages -d
chmod 600 id_rsa_gh-pages
eval `ssh-agent -s`
ssh-add id_rsa_gh-pages

# Now that we're all set up, we can push.
git push $SSH_REPO $TARGET_BRANCH

rm id_rsa_gh-pages

