import re
import sys
import fileinput
import os

# Set this to the name of the repo, if you don't want it to be read from the filesystem.
# It assumes the changelog file is in the root of the repo.
repo_name = ""

# This script goes through the provided file, and replaces any " \#<number>",
# with the valid mark down formatted link to it. e.g.
# " [\#number](https://github.com/arkworks-rs/template/pull/<number>)
# Note that if the number is for a an issue, github will auto-redirect you when you click the link.
# It is safe to run the script multiple times in succession. 
#
# Example usage $ python3 linkify_changelog.py ../CHANGELOG.md
if len(sys.argv) < 2:
    print("Must include path to changelog as the first argument to the script")
    print("Example Usage: python3 linkify_changelog.py ../CHANGELOG.md")
    exit()

changelog_path = sys.argv[1]
if repo_name == "":
    path = os.path.abspath(changelog_path)
    components = path.split(os.path.sep)
    repo_name = components[-2]

for line in fileinput.input(inplace=True):
    line = re.sub(r"\- #([0-9]*)", r"- [\\#\1](https://github.com/arkworks-rs/" + repo_name + r"/pull/\1)", line.rstrip())
    # edits the current file
    print(line)