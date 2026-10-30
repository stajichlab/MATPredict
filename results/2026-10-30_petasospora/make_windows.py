"""Extract the four MAT windows (DIC1 start -2 kb to SLA2 end +1 kb) into data/windows.fa; header: name contig start end."""
import subprocess, sys
D = sys.argv[1]  # directory with the four .fa
W = [("GCA_022641835.1_ASM2264183v1", "CP048715.1", 626500, 645500, "pf_a_cp048715"),
     ("GCA_001599195.1_JCM_3601_assembly_v001", "BCGI01000001.1", 639000, 658500, "pf_a_jcm3601"),
     ("GCA_001983305.1_ASM198330v1", "MPUK01000005.1", 609500, 631500, "pf_alpha_strain65"),
     ("GCA_050575385.1_ASM5057538v1", "CM115473.1", 629500, 652000, "pf_alpha_tbrc4498")]
with open("data/windows.fa", "w") as o:
    for g, c, a, b, n in W:
        s = "".join(subprocess.check_output(["samtools", "faidx", f"{D}/{g}.fa", f"{c}:{a}-{b}"], text=True).split("\n")[1:])
        o.write(f">{n} {c} {a} {b}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
