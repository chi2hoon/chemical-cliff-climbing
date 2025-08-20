few_shot_examples = [
    {
        "input": {
            "mol1": "c1ccccc1",
            "mol2": "Cc1ccccc1",
            "diff_type": "lipophilicity"
        },
        "output": (
            "Major Structural Change: Toluene replaces benzene's hydrogen with a methyl group.\n"
            "Physicochemical Impact: The methyl increases hydrophobic surface and lipophilicity, slightly raising steric bulk.\n"
            "Biological Activity Impact: Enhanced hydrophobic contacts can strengthen binding in a nonpolar pocket, potentially improving potency; sterics may also bias preferred conformations.\n"
            "Key Interactions: Increased van der Waals contacts from the methyl substituent within hydrophobic regions."
        )
    },
    {
        "input": {
            "mol1": "Cc1ccccc1",
            "mol2": "Cc1ccncc1",
            "diff_type": "electronic effect"
        },
        "output": (
            "Major Structural Change: Ring nitrogen introduced (phenyl → pyridyl).\n"
            "Physicochemical Impact: Increased polarity and H-bond basicity; altered electron distribution reduces ring hydrophobicity.\n"
            "Biological Activity Impact: Possible new H-bond acceptor interactions; however, reduced hydrophobicity may weaken binding in lipophilic sites.\n"
            "Key Interactions: Potential H-bond to protein donor near ring nitrogen; altered π-stacking characteristics."
        )
    },
    {
        "input": {
            "mol1": "CC(C)O",
            "mol2": "CC(C)OC",
            "diff_type": "steric/electronic"
        },
        "output": (
            "Major Structural Change: Alcohol → methyl ether substitution.\n"
            "Physicochemical Impact: Reduced H-bond donor capacity; slight increase in lipophilicity and steric bulk.\n"
            "Biological Activity Impact: Loss of donor interactions may reduce affinity if a donor is required; increased lipophilicity can improve membrane permeation.\n"
            "Key Interactions: Loss of direct H-bond donation; potential maintenance of acceptor capability via ether oxygen."
        )
    }
]
