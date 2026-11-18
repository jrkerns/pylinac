# Divergence from jrkerns/pylinac
Based on upstream 3.50.0 (5c6a17c0). Last synced: 2026-11-18.

## Upstreamable
| Change | Branch/commits | Status |
|---|---|---|
| Cherry pick changes from DCM-stackorder-551 | 551-DCM-stackorder-revised |  |
| Flatness of 0 treated as int | fix-flat-int | not yet submitted |

## Fork-specific (do not upstream)
| Change | Reason |
|---|---|
| Removed GUI dependencies | LinaQA bundles its own GUI |
| Changed results() return format | Needed by LinaQA, breaks API |
