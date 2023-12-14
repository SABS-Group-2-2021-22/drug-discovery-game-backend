import requests
from requests.exceptions import HTTPError

LLM_API_URL = "http://127.0.0.1:8080/completion"

def query(payload):
    response = requests.post(LLM_API_URL, json=payload)
    return response.json()

def process_prompt(prompt):
    llm_input = {
        "prompt": prompt,
        "n_predict": 200
    }
    response = query(llm_input)
    print(response)
    return response["content"]

