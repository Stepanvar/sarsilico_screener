from django.core.management.base import BaseCommand
from app.models import PredefinedProteinTarget

class Command(BaseCommand):
    help = 'Populate the database with predefined protein targets'

    def handle(self, *args, **kwargs):
        targets_data = [
            {
                "name": "Nonstructural protein 1 (nsp1)",
                "description": (
                    "Nsp1 is shown to promote cellular mRNA degradation, block host cell translation, "
                    "and inhibit the innate immune response to virus infection. The crystal structure of "
                    "nsp1 globular domain (residue 13-121) was downloaded from the PDB database with code 7K7P "
                    "and prepared for peptide/antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 2 (nsp2)",
                "description": (
                    "Nsp2 is a protein containing three zinc fingers, indicating a RNA binding site. The full length "
                    "structure, suggesting its role in linking viral transcription within the RTC to the translation "
                    "initiation of the viral message, was obtained from 7MSW (residues 1-638) and from 7EXM (residues 1-277) "
                    "for the N-terminal zinc ion binding domain, provided for peptide or antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 3 (nsp3)",
                "description": (
                    "Nsp3 is a large multi-domain protein. Its ADP-ribose phosphatase domain (ADRP/MacroD) is involved "
                    "in the host immune response. The structure of nsp3 MacroD in complex with AMP and MES was downloaded from "
                    "PDB with code 6W6Y. Two binding sites (AMP_site and MES_site) are defined. Additional structures from 5RSF, "
                    "7KAG, 7LGO, and 7RQG are provided for docking."
                )
            },
            {
                "name": "Papan-like protease (PLpro)",
                "description": (
                    "PLpro cleaves the nsp1/2, nsp2/3, and nsp3/4 boundaries and works with Mpro to process polyproteins. "
                    "Crystal structures of wild type and C111S mutant in complex with a compound (codes 7RZC and 7SQE) are provided "
                    "for small molecule and peptide/antibody docking."
                )
            },
            {
                "name": "Main protease (Mpro, nsp5)",
                "description": (
                    "Also known as chymotrypsin-like protease (3CLpro), Mpro cleaves the majority of sites in polyproteins. "
                    "The room-temperature X-ray structure of Mpro in complex with PF-07321332 (code 7SI9) is prepared for small molecule "
                    "and peptide/antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 6 (nsp6)",
                "description": (
                    "No experimental structure is currently available. A computationally predicted structure from Zhang’s lab was "
                    "used and prepared for peptide or antibody docking."
                )
            },
            {
                "name": "Nonstructural protein12/7/8 (nsp12/7/8, RdRp)",
                "description": (
                    "Nsp12 is the RNA-dependent RNA polymerase that complexes with nsp7 and nsp8. Structures include a complex with "
                    "RNA and Remdesivir (code 7BV2) defining the RTP binding site, and a structure with suramin (code 7D4F) for the RNA binding site. "
                    "For peptide/antibody docking, both the complex and individual chains (nsp12, nsp7, and nsp8) are provided."
                )
            },
            {
                "name": "Nonstructural protein 9 (nsp9)",
                "description": (
                    "Nsp9 may function as a single-stranded RNA binding protein. Crystal structures in dimer (6WXD) and monomer (6W9Q) forms "
                    "are provided, along with small molecule docking files targeting the oligomerization interface (7KRI) and a conserved site near "
                    "the C-terminal GxxxG-helix (7N3K)."
                )
            },
            {
                "name": "Nonstructural protein 10 (nsp10)",
                "description": (
                    "Nsp10 stimulates the exoribonuclease and 2'-O-methyltransferase activities of nsps 14 and 16 by forming complexes with them. "
                    "Structures include nsp10/14 and nsp10/16 complexes (codes 7ORR, 7ORU) and the unbound form (6ZPE) for peptide/antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 10/14 (nsp10/14)",
                "description": (
                    "Nsp14 is bifunctional, containing an N-terminal exoribonuclease (ExoN) domain and a C-terminal N7-methyltransferase (N7-MTase) domain. "
                    "Stabilized by nsp10, its docking files for the ExoN and N7-MTase sites (and a chapso binding site) are prepared from the Cryo-EM structure (7N0D). "
                    "Both the complex and the extracted nsp14 structure are provided for peptide/antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 16/10 (nsp16/10)",
                "description": (
                    "Nsp16 is a SAM-dependent 2'-O-methyltransferase that is active only when bound to nsp10. Structures include the complex with "
                    "7-methyl-GpppA, SAM, and MGP (code 6WVN) and a higher resolution complex with SAM (code 6W4H). Docking files are provided for SAM_site, "
                    "GTA_site, and MGP_site, as well as for peptide/antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 13 (nsp13, helicase)",
                "description": (
                    "The helicase unwinds duplex oligonucleotides in an NTP-dependent manner. Its structure in complex with an ATP analog (7NN0) is provided "
                    "for small molecule docking at the ANP site, along with additional docking files for a fragment binding site (5RML) and both apo (7NIO) and "
                    "bound forms for peptide/antibody docking."
                )
            },
            {
                "name": "Nonstructural protein 15 (nsp15)",
                "description": (
                    "Nsp15 is a uridylate-specific endoribonuclease implicated in evading the innate immune response. Structures include a complex with "
                    "Tipiracil (7K1L) for small molecule docking and both monomer (6VWW) and hexamer (7N06) forms for peptide/antibody docking."
                )
            },
            {
                "name": "ORF 3A",
                "description": (
                    "Accessory protein ORF3A dimer structure (code 7KJR) is provided for peptide and antibody docking."
                )
            },
            {
                "name": "ORF 7A",
                "description": (
                    "Accessory protein ORF7A crystal structure (code 7CI3) is provided for peptide or antibody docking."
                )
            },
            {
                "name": "ORF 8",
                "description": (
                    "Accessory protein ORF8 crystal structure (code 7JX6) is provided for peptide and antibody docking."
                )
            },
            {
                "name": "ORF 9B",
                "description": (
                    "Accessory protein ORF9B crystal structure (code 6Z4U) is provided for peptide and antibody docking."
                )
            },
            {
                "name": "Spike protein (S protein)",
                "description": (
                    "The spike protein is a trimer composed of S1-S2 heterodimers. The receptor binding domain (RBD) from the complex with ACE2 "
                    "(code 6M0J) is provided for peptide/antibody docking. Additionally, the open (6VYB) and closed (6ZGE) states of the trimer are available."
                )
            },
            {
                "name": "S2 of S protein",
                "description": (
                    "The S2 segment of the spike protein forms a post-fusion 6-helical bundle (code 7COT), provided for peptide or antibody docking."
                )
            },
            {
                "name": "N-terminal domain (NTD) of S protein",
                "description": (
                    "The NTD, a dominant epitope for antibody binding, is provided from both the wild type (7B62) and Kappa variant (7SOD) structures."
                )
            },
            {
                "name": "Envelop small membrane protein (E protein)",
                "description": (
                    "The envelope protein forms a pentameric ion channel (E channel). Its modeled structure is provided for peptide or antibody docking."
                )
            },
            {
                "name": "Membrane protein (M protein)",
                "description": (
                    "The membrane protein is involved in coronavirus assembly and acts as a protective antigen. The structure (code 8CTK) is provided "
                    "for peptide or antibody docking."
                )
            },
            {
                "name": "Nucleocapsid protein (N protein)",
                "description": (
                    "The nucleocapsid protein, which forms a ribonucleoprotein complex with viral RNA, is provided in multiple forms: N-terminal domain (6YI3), "
                    "C-terminal domain (6YUN), and full-length (8FD5) for peptide/antibody docking. A ribonucleotide-binding site (NCB site) is also prepared."
                )
            },
            {
                "name": "Angiotensin-converting enzyme 2 (ACE2)",
                "description": (
                    "The human ACE2 receptor, extracted from the spike RBD complex (6M0J), is provided for peptide or antibody docking."
                )
            },
        ]

        for data in targets_data:
            target, created = PredefinedProteinTarget.objects.get_or_create(
                name=data["name"],
                defaults={"description": data["description"]}
            )
            if created:
                self.stdout.write(f"Created target: {target.name}")
            else:
                self.stdout.write(f"Target already exists: {target.name}")
