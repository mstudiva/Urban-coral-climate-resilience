import re,csv,collections
H="Heat shock proteins"; O="Oxidative stress/antioxidants"; A="Apoptosis/cell death"; D="DNA repair"; I="Immune response"; U="Ubiquitin-proteasome"; S="Other stress-related"
# ordered rules: first match wins. Patterns are matched case-insensitively against the annotation text.
RULES=[
 # --- exclusions (no family), checked first: ubiquitin-like modifiers other than ubiquitin, generic domains, vertebrate
 #     lymphocyte/cytokine GO terms and long descriptions that only mention a keyword in passing
 (None, r"ubiquitin-like protein-specific protease|SUMO|sumoylation|small ubiquitin|neddyl|NEDD8|ISG15|ubiquitin-like modifier|ubiquitin homologues|Rad60"
        r"|histone acetylation \(HAT\) complex SAGA|Involved in intracellular signal transduction mediated by cytokines|Required for RNA-mediated gene silencing"
        r"|Transcriptional regulator\. Recognizes and binds to the DNA sequence|brca1 gene 1|Component of the BLOC-1 complex|Concanavalin A-like"
        r"|B cell (affinity|differentiation)|germinal center|T cell (differentiation|polarity)|CD8-positive|establishment of T cell|Immunoglobulin (C-2|domain)|immunoglobulin domain"
        r"|Interleukin enhancer-binding|Interleukin-like EMT|interleukin-\d+ (production|biosynthetic)|regulation of interleukin|Interleukin 2 receptor"
        r"|Leukocyte cysteine proteinase|macrophage erythroblast attacher|fat cell proliferation|into host cell cytoplasm|iron-sulfur cluster co-chaperone|cytochrome P450"),
 # --- immune signalling that would otherwise match the apoptosis rule
 (I, r"TNF receptor-associated factor|\bTRAF\d?\b|lipopolysaccharide-induced|LPS-induced tumor necrosis factor|Tumor necrosis factor, alpha-induced protein 3"),
 # --- oxidative stress that would otherwise match another rule
 (O, r"toxic effects of hydrogen peroxide"),
 # --- heat shock proteins and their co-chaperones
 (H, r"heat shock|\bHSP ?\d|hsp20|hsp40|hsp70|hsp90|hsc70|\bGrp94\b|DnaJ|\bDnaK\b|chaperonin|BAG family molecular chaperone|suppressor of G2 allele of SKP1|Hsp70 interacting|Hsp70 protein binding|alpha[- ]crystallin|T-complex protein 1|assists the folding of proteins upon ATP hydrolysis|stress-induced-phosphoprotein"),
 # --- DNA repair
 (D, r"DNA repair|base-excision repair|nucleotide-excision repair|mismatch repair|double-strand break repair|interstrand cross-link repair|DNA damage|damaged DNA|\bPMS1\b|\bMSH\d|\bMLH\d|mutS protein|DNA polymerase zeta|uracil[- ]DNA|Excises uracil|DNA glycosylase|\bRAD\d\d|\bXRCC|\bERCC|\bXP-?[ACG]\b|XPA binding|photo-?lyase|translesion|Fanconi|ATP-dependent DNA helicase|RecQ|UvrD|Helicase-like transcription factor|DNA annealing helicase|flap endonuclease|poly ?\(ADP-ribose\) polymerase|\bPARP|BRCA1-A complex|AP endonuclease|deoxyribodipyrimidine|alkylated DNA|methylated-DNA|8-oxo|non-homologous end joining"),
 # --- apoptosis / cell death
 (A, r"apopto|caspase|\bBcl-?2|\bBIR domain|\bBIRC|inhibitor of apoptosis|death domain|death effector|death-associated|death receptor|programmed cell death|tumou?r necrosis factor|\bTNF\b|\bTNFR|MDM2|\bp53\b|necroptosis|pyroptosis|phosphatidylserine exposure"),
 # --- ubiquitin-proteasome
 (U, r"ubiquitin|proteasom|RING finger and|ring finger and CHY|\bcullin|\bF-box|SCF-dependent|\bUBX domain|deubiquitin|Kelch-like protein 12\b|HECT-domain|UBR box|LON peptidase N-terminal domain and RING"),
 # --- oxidative stress / antioxidants
 (O, r"glutathione|glutamate-cysteine ligase|superoxide|\bcatalase\b|peroxiredoxin|thioredoxin|glutaredoxin|peroxidase|peroxidasin|antioxidant|oxidative stress|Oxidation resistance|redox homeostasis|reactive oxygen species|\bferritin\b|Stores iron in a soluble|nitric[- ]oxide synthase|nitric oxide (catabolic|biosynthetic)|methionine \(S\)-S-oxide reductase|methionine sulfoxide reductase|sulfiredoxin|protein[- ]disulfide (reductase|oxidoreductase)|selenoprotein|amine oxidase"),
 # --- immune response
 (I, r"immun|complement (activation|component|factor|control)|SUSHI repeat|properdin|anaphylatoxin|MAC/Perforin|lectin|toll[- ]like|\bTIR domain|interferon|NACHT|NF-kappa|Nuclear factor of kappa|scavenger receptor|tachylectin|mannose receptor|fibrinogen-related|defense response|antimicrobial|bactericidal|lysozyme|MyD88|tyrosinase|AIG1 family|macrophage|leukocyte migration|T cell proliferation involved in immune|chemokine|cytokine"),
 # --- other stress-related
 (S, r"universal stress protein|stress response|response to stress|stress-induced|stress-activated|unfolded protein|ER overload|endoplasmic reticulum stress|ER-associated misfolded|mis-folding|eukaryotic translation initiation factor 2-alpha kinase|hypoxia|Egl-9|MAP kinase|mitogen-activated protein kinase|peptidyl-prolyl cis-trans isomerase|PPIases accelerate|\bFKBP|protein disulfide isomerase|protein folding|calreticulin|calnexin|prefoldin|plasma membrane repair|osmotic stress|autophag"),
]
RULES=[(f,re.compile(p,re.I)) for f,p in RULES]
def classify(annot):
    if annot is None or annot in ("","NA","-"): return None
    for fam,rx in RULES:
        if rx.search(annot): return fam
    return None
