# Contributing

Thank you for considering making contributions to `arkworks-rs/groth16`!

Contributing to this repo can be done in several forms, such as participating in discussion or proposing code changes. 
To ensure a smooth workflow for all contributors, the following general procedure for contributing has been established:

1) Either open or find an issue you'd like to help with
2) Participate in thoughtful discussion on that issue
3) If you would like to contribute:
    * If the issue is a feature proposal, ensure that the proposal has been accepted
    * Ensure that nobody else has already begun working on this issue. 
    If they have, please try to contact them to collaborate
    * If nobody has been assigned for the issue and you would like to work on it, make a comment on the issue to inform the community of your intentions to begin work. (So we can avoid duplication of efforts)
    * We suggest using standard Github best practices for contributing: fork the repo, branch from the HEAD of master, make some commits on your branch, and submit a PR from the branch to master.
    More detail on this is below
    * Be sure to include a relevant change log entry in the Pending section of CHANGELOG.md (see file for log format)
        * If the change is breaking, we may add migration instructions.

Note that for very small or clear problems (such as typos), or well isolated improvements, it is not required to an open issue to submit a PR.
But be aware that for more complex problems/features touching multiple parts of the codebase, if a PR is opened before an adequate design discussion has taken place in a github issue, that PR runs a larger likelihood of being rejected.

Looking for a good place to start contributing? How about checking out some good first issues

## Branch Structure

`groth16` has its default branch as `master`, which is where PRs are merged into. Releases will be periodically made, on no set schedule.
All other branches should be assumed to be miscellaneous feature development branches.

All downstream users of the library should be using tagged versions of the library pulled from cargo.

## How to work on a fork
Please skip this section if you're familiar with contributing to opensource github projects.

First fork the repo from the github UI, and clone it locally.
Then in the repo, you want to add the repo you forked from as a new remote. You do this as:
```bash
git remote add upstream git@github.com:arkworks-rs/groth16.git
```

Then the way you make code contributions is to first think of a branch name that describes your change.
Then do the following:
```bash
git checkout master
git pull upstream master
git checkout -b $BRANCH_NAME
```
and then work as normal on that branch, and pull request to upstream master when you're done =)

## Updating documentation

All PRs should aim to leave the code more documented than it started with.
Please don't assume that its easy to infer what the code is doing, 
as that is usually not the case for these complex protocols. 
(Even when you understand the paper!)

Its often very useful to describe what is the high level view of what a code block is doing,
and either refer to the relevant section of a paper or include a short proof/argument for why it makes sense before the actual logic.

## Performance improvements

All performance improvements should be accompanied with benchmarks improving, or otherwise have it be clear that things have improved.
For some areas of the codebase, performance roughly follows the number of field multiplications, but there are also many areas where
hard to predict low level system effects such as cache locality and superscalar operations become important for performance. 
Thus performance can often become very non-intuitive / diverge from minimizing the number of arithmetic operations.