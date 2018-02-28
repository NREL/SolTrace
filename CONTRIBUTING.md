# Contributing to SolTrace

You can contribute to SolTrace by letting us know about problems or suggesting new features, or by making your own changes or additions to the code. You may want to help us fix an [issue someone else reported](https://github.com/NREL/SolTrace/issues), fix an issue you discovered, or add a new feature to SolTrace.

## Let us know about a problem or suggest a new feature

If you find a problem with SolTrace, or would like to request a new feature, let us know by [submitting an issue](https://github.com/NREL/SolTrace/issues/new).

If you have a question about using SolTrace, you can ask us at [SolTrace support](mailto://soltrace.support@nrel.gov).

## Contribute code to the SolTrace project

If you are ready to contribute code to SolTrace, there are a couple of things you should know first:

* First off, please read [SolTrace's Contribution Policy](https://github.com/NREL/SolTrace/wiki/Contribution-Policy).  In particular, you'll need to send us an email (to [soltrace.support@nrel.gov](mailto://soltrace.support@nrel.gov)) that states:

_I agree to contribute to SolTrace. I agree to the following terms and conditions for my contributions: First, I agree that I am licensing my contributions under the terms of the current SolTrace license. Second, I agree that, in order to conform to any future open source software license(s) under which SolTrace may be provided, the terms of my license may be modified without any notice to me and without my consent. Third, I represent and warrant that I am authorized to make the contributions and grant the license. If my employer has rights to intellectual property that includes my contributions, I represent and warrant that I have received permission to make contributions and grant the required license on behalf of my employer._

* We use GitHub to manage the open source project, so you will need to learn how to use it to fork, clone, branch, check out, pull, add, commit, and push your work. 

### Instructions

Here are the steps we would like you to follow when you contribute code to SAM:

1. Install GitHub on your computer.
2. Follow the instructions on the [SolTrace wiki](https://github.com/NREL/SolTrace/wiki) to clone the SolTrace repositories and build SolTrace.
3. Create a fork on GitHub.com for the repository (SolTrace, SolTrace, SSC, LK, or WEX) you are contributing to.
4. Clone your fork to create a local copy on your computer.
5. Create a branch for your changes.
6. Make your changes to the code.
7. Build SolTrace and test it to make sure your code works as expected (see [below](#test-protocol)).
8. Commit and push the changes to the branch.
9. Create a pull request for us to review the branch. If the changes meet our requirements, we will merge the branch into the main repository.

### Resources for Learning GitHub

If you are new to GitHub, you can find helpful articles to help you learn how it works on the web. Some examples are:

* [Using the Fork-and-Branch Git Workflow](https://blog.scottlowe.org/2015/01/27/using-fork-branch-git-workflow/) from Scott's Weblog is a good introduction.

* [Git Concepts: Branches, Forks, and Pull Requests](http://willi.am/blog/2014/05/12/git-concepts-branches-forks-and-pull-requests/) from Will Anderson is useful, although the video on the page does not work.

* [3.2 Git Branching - Basic Branching and Merging](https://www.git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) from the Git documentation.

* [Fork a Repo](https://help.github.com/articles/fork-a-repo/) from GitHub Help.

* [About pull requests](https://help.github.com/articles/about-pull-requests/) from GitHub Help.

### Test Protocol

We are in the process of setting up a Google Test framework for testing your contribution to ensure that it does not cause any problems with the software. 

For now, you can help to ensure that your code works with the rest of SolTrace by:

1. Compiling SolTrace with your contribution for Windows, Mac, and Linux.

2. Fixing any compiler warning messages.

3. Running the compiled program with several sample files.

**Note that a change to _coretrace_ may affect performance in [SolarPILOT](https://github.com/NREL/SolarPILOT). Changes to the ray tracing code (not the interface) should be tested within both SolTrace and SolarPILOT.** 
