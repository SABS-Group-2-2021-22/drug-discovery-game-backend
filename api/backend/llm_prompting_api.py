import requests
from requests.exceptions import HTTPError


def process_prompt(prompt):

    LLM_API_URL = "https://api-inference.huggingface.co/models/meta-llama/Llama-2-13b-chat-hf"
    api_key = "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"
    
    headers = {"Authorization": f"Bearer {api_key}"}
    def query(payload):
        response = requests.post(LLM_API_URL, headers=headers, json=payload)
        return response.json() 
    response = query({
        "inputs": prompt,
    })
    print(response)
    return response

