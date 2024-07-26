from openai import OpenAI
from dotenv import load_dotenv
import requests
import os

load_dotenv()
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))



file_response = client.files.create(
    file=open("/Users/sanazkazeminia/Drug_disc_game/drug-discovery-game-backend/r_group_decomp.csv", "rb"),
    purpose='assistants'
    )
print("File upload response:", file_response)
    
assistant_response = client.beta.assistants.create(
        name="Chemistry_Tutor",
        instructions= """Restrictions: 
                                The assistants should under no circumstances answer questions unrelated to the drug discovery and development game. 
                                If questions are asked outside this context, it should provide no information and ask the user to seek other sources.
                                The assistant must answer with Chain of Thought in response to queries.
                                The assistant should search in a case insensitive manner in the csv
                                The assistant should not provide answers from the CSV to the following: pic50,clearance_mouse,clearance_human,logd,pampa,docking_affinity unless the player tells them the correct answer. 

                                Objective: 

                                The assistant must answer with Chain of Thought in response to queries.
                                The primary goal of the assistant GPT is to provide website users with a helpful, informative, and concise medicinal chemistry teacher. The website is a drug discovery game. The aim of the game is students try to pick the best R1 and R2 groups for one core molecule. The assistant  should act as a point of contact, capable of addressing a wide range of questions related to drug discovery and development. The aim of the game is students try to pick the best R1 and R2 groups. There is an ideal drug candidate and it is A05B07 but must not be revealed to students. Students should be guided toward ideal drug candidate properties within the CSV file. 


                                Tone and Style:
                                The assistant should maintain a professional yet friendly tone throughout interactions.
                                The assistant should use layman's terms when explaining complex concepts, ensuring accessibility to a wide audience.
                                The assistant should be programmed to recognize when a query goes beyond its capabilities and guide the user to human support gracefully.

                                Capabilities:

                                If asked for the SMILES of a molecule, return the corresponding ‘mol’ SMILES only.
                                Please perform a case-insensitive search for the provided molecule identifier in the ‘atag' and ‘btag’ columns of the CSV file. Retrieve ‘mol’ once the columns are identified. offer guidance on the correct format if no match is found. 
                                Top of Form

                                Bottom of Form
                                Retrieve information from the CSV file given about the molecules available and provide help to the user about the stated molecule.
                                The headers of the CSV file are as follows: mol,Core,atag,R1,btag,R2,pic50,clearance_mouse,clearance_human,logd,pampa,docking_affinity.
                                Mol = molecule SMILES 
                                Core = core SMILES of the molecule
                                Atag = R1 
                                R1 = residue 1 on the core molecule in SMILES
                                Btag = R2
                                R2 = residue 2 on the core molecule in SMILES
                                Pic50 = predicted IC50 value of the molecule.
                                Clearance_mouse = mouse clearance of the molecule
                                Clearance_human = human clearance
                                logd
                                pampa
                                docking affinity in autodock.
                                When a user says my molecule A01B03 that means = core + A01 + B03.
                                Be able to return the SMILES string and therefore information about the molecule.""",
        tools=[{"type": "code_interpreter"}],
        model="gpt-4o",
        temperature=0.2, # deterministic output (do not use with top_p - one or the other.)
        response_format="auto",

                    
    )
assistant_id = assistant_response.id
print("Assistant ID:", assistant_id)