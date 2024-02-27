from openai import OpenAI
from dotenv import load_dotenv
import requests
load_dotenv()

client = OpenAI()

try:
    file_response = client.files.create(
      file=open("/Users/sanazkazeminia/Drug_disc_game/drug-discovery-game-backend/r_group_decomp.csv", "rb"),
      purpose='assistants'
    )
    print("File upload response:", file_response)
    
    assistant_response = client.beta.assistants.create(
        name="GAME_ASSISTANT",
        instructions= "The primary goal of the assistant GPT is to provide website users with a helpful, informative, and concise medicinal chemistry teacher. The website is a drug discovery game. The aim of the game is students try to pick the best R1 and R2 groups for one core molecule. The assistant  should act as a point of contact, capable of addressing a wide range of questions related to drug discovery and development. The aim of the game is students try to pick the best R1 and R2 groups. There is an ideal drug candidate and it is A05B07 but it must not be revealed to students. Students should be guided toward ideal drug candidate properties within the CSV file.",
        tools=[{"type": "code_interpreter"}],
        model="gpt-4-turbo-preview"
    )
    
    thread_response = client.beta.threads.create()
    print("Thread creation response:", thread_response)
    
    message_response = client.beta.threads.messages.create(
        thread_id=thread_response.id,
        role="user",
        content="what is logD? how do improve it?"
    )    
    run_response = client.beta.threads.runs.create(
      thread_id=thread_response.id,
      assistant_id=assistant_response.id,
      instructions="Please address the user as Jane Doe. The user has a premium account."
    )
    
    messages_response = client.beta.threads.messages.list(
      thread_id=thread_response.id
    )
    print("Messages list response:", messages_response)
except Exception as e:
    print("An error occurred:", e)

import requests

# Replace these variables with your actual thread ID and run ID
thread_id = "thread_lOM7RvF33IgeEJUyreO9igRY"
run_id = "run_dIHtncewbt400g2ojRA3nGyw"

url = f"https://api.openai.com/v1/threads/{thread_id}/runs/{run_id}"

headers = {
    "Authorization": f"Bearer {}"
}

response = requests.get(url, headers=headers)

if response.status_code == 200:
    run_details = response.json()
    print("Run Details:", run_details)
    # Extract and print the status
    status = run_details.get("status", "No status found")
    print("Run Status:", status)
else:
    print(f"Failed to retrieve run details. Status code: {response.status_code}")

