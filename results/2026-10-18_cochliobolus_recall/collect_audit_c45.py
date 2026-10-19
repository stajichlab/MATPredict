import shutil, os
S = '/private/tmp/claude-501/-Users-jstajich-projects-MATPredict/225e2914-af32-4993-9a12-4f0df8da228b/scratchpad/diag'
D = 'results/2026-10-18_cochliobolus_recall/audit_run_c45'
for g in os.listdir(S):
    p = os.path.join(S, g)
    if os.path.isdir(p):
        os.makedirs(f'{D}/{g}', exist_ok=True)
        for f in ('detection_report.yaml', 'evidence_diagnostics.jsonl'):
            shutil.copy(f'{p}/{f}', f'{D}/{g}/{f}')
