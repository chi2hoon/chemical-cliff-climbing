# File: activity_cliff_app/prompts.py

import os
from openai import OpenAI
from rag.retriever import retrieve_examples
from rdkit import Chem
from structure_diff import detect_diff_type
from prompt_examples import few_shot_examples

SYSTEM_PROMPT = """You are an expert computational chemist in the field of drug discovery. You must clearly and concisely explain the reasons for the structural differences between two compounds and the resulting changes in their biological activity. Your explanation should follow this structure:

**Major Structural Change:**
[Describe the main structural changes between the two compounds. E.g., substitution, removal, stereochemical changes, etc.]

**Physicochemical Impact:**
[Explain the impact of these structural changes on the molecule's physicochemical properties (e.g., hydrophobicity, polarity, electron density, etc.).]

**Biological Activity Impact:**
[Explain the potential impact of physicochemical property changes on biological activity (e.g., target binding, metabolic stability, absorption, etc.).]

**Key Interactions:**
[If applicable, mention specific interactions (e.g., hydrogen bonding, hydrophobic interactions, steric hindrance) that help explain the activity change. If not applicable, omit this section or state 'Not Applicable'.]

The explanation should be written in professional and easy-to-understand language. Focus on core information without unnecessary introductions or conclusions.
"""

def generate_prompt(row):
    user_prompt = f"""Two compounds:
- SMILES 1: {row['mol1_smiles']} (pIC50: {row['mol1_activity']:.2f})
- SMILES 2: {row['mol2_smiles']} (pIC50: {row['mol2_activity']:.2f})

Explain the reasons for the activity change due to structural differences between the above two compounds, following the given guidelines."""
    return SYSTEM_PROMPT, user_prompt

def generate_fewshot_prompt(row):
    shots = ""
    for ex in few_shot_examples:
        # Assuming few_shot_examples already contain the '자동화된 해석 및 가설' format
        shots += f"Q: {ex['input']['mol1']} vs {ex['input']['mol2']}\nA: {ex['output']}\n\n"
    user_prompt = shots + generate_prompt(row)[1] # Append to the user part of the basic prompt
    return SYSTEM_PROMPT, user_prompt

def generate_rag_prompt(row, client):
    mol1 = Chem.MolFromSmiles(row["mol1_smiles"])
    mol2 = Chem.MolFromSmiles(row["mol2_smiles"])
    diff_type = detect_diff_type(mol1, mol2)
    examples = retrieve_examples(client, row["mol1_smiles"], row["mol2_smiles"], diff_type)
    body = "\n".join([f"- 예시: {ex}" for ex in examples])

    rag_context = f"The following are examples of similar structural and activity changes. Refer to this information to formulate your answer:\n{body}\n\n"

    user_prompt = f"""Two compounds:
- SMILES 1: {row['mol1_smiles']} (pIC50: {row['mol1_activity']:.2f})
- SMILES 2: {row['mol2_smiles']} (pIC50: {row['mol2_activity']:.2f})

Structural Difference Type: {diff_type}

Explain the reasons for the activity change due to structural differences between the above two compounds, referring to the given guidelines and provided similar examples.
"""
    return SYSTEM_PROMPT, rag_context + user_prompt

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

def call_llm(system_prompt, user_prompt):
    response = client.chat.completions.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        temperature=0.3
    )
    return response.choices[0].message.content.strip()

def translate_text(text, target_language="Korean"):
    translation_client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
    response = translation_client.chat.completions.create(
        model="gpt-4o-mini", # Using a cheaper model for translation
        messages=[
            {"role": "system", "content": f"You are a professional translator. Translate the following text into {target_language}. Ensure technical terms are translated accurately and contextually."},
            {"role": "user", "content": text}
        ],
        temperature=0.1
    )
    return response.choices[0].message.content.strip()