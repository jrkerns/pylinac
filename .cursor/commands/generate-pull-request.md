# generate-pull-request

Generating a Pull Request is against the pylinac BitBucket repository.
PRs should be created in "draft" mode, not in "ready" mode.
Unless otherwise stated, the PR should be against the "master" branch.

## Git

The branch name should ideally include the JIRA ticket, but is not required. Prompt me
if it does not contain a ticket reference before creating a PR; it's possible there isn't one.

Examples:
* docs/RAM-1234-add-reference-to-requirements
* feature/RAM-2345-add-standard-profile-format
* bugfix/RAM-3456-fix-description-injection-issue

## BitBucket

### PR Titles

PR titles should be the ticket and a short summary, similar to the branch name.
Examples:
* RAM-4567 Proof of concept of using Python .whl files
* RAM-5678 Use Docker:28, GCP can't handle Docker:29

### PR Descriptions

PR descriptions should include a ticket and Summary. Screenshots are also a nice-to-have for user-facing features.

Example 1:

Ticket
RAM-1144

Summary

This is the SETUP for fixing cybersecurity issues w/ the sandbox. This adds pen test scripts to calculation tests. Currently, they are set to wrong values, but the intent is to do a multi-layered approach to the fixes, using this PR/branch as a baseline. I.e. fork from this branch with strategy 1, ensure it works, fork from this w/ strategy 2, ensure it works, etc. This proves that we currently have a security hole and reproduces it.

Example 2:

Ticket
n/a

Summary

Follow up PR to the parser backend to enable it on the front end.

Small todo left is to only add xlsx to the supported file types if the profile data flag is enabled. I based this PR on master which doesn’t have the flag yet, but I’ll update the PR after the UI prs are merged.

### PR Reviewers

There should be 2 reviewers for PRs.

### Comments

Comments by the PR author are highly encouraged. These are done to explain
logic that may otherwise be confusing. Unlike comments within the source code,
comments on the code within the BitBucket Pull request may be relevant only for the
PR review. For example, if we switch from one methodology to another, that
is something that would make for an excellent comment at line(s) where the change is at.
