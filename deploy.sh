#!/bin/bash
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
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH

# Clean out existing contents
#rm -rf out/**/* || exit 1
find . -mindepth 1 -maxdepth 1 ! -name '.git' ! -name "README" | xargs rm -rf 

# Now let's go have some fun with the cloned repo
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"

#Main code

#Download index.html
wget -nv https://github.com/orex/supercell/raw/deploy/index.html

mkdir -p doc
mkdir -p exe

cp ${c_path}/build/doc/man/supercell_man.pdf doc/.
cp ${c_path}/build/doc/man/supercell_man.html doc/.
cp ${c_path}/doc/man/supercell_man.css doc/.
cp ${c_path}/build/src/sc_cli/supercell exe/.
cd exe
tar czvf supercell.tar.gz supercell
rm supercell

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add --no-ignore-removal .
git commit -m "Deploy to GitHub Pages: ${SHA}"

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

