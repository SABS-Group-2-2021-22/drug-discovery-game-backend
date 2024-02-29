from openai import OpenAI
import time
from dotenv import load_dotenv
import os


load_dotenv()
client = OpenAI()

assistant_id = os.getenv("ASSISTANT_ID")

try:
    #thread_response = client.beta.threads.create() # this creates a thread for which the ID can be used to continue the conversation 
    #print("Thread creation response:", thread_response)
        
    message_response = client.beta.threads.messages.create(
            thread_id='thread_Jj24BQoMxDGWvoZw4F7LoeFj',
            role="user",
            content="what is the molecule SMILES for A17B12 and how would i improve its pampa?"
        )    
    print("Message creation response:", message_response)
    
    run_response = client.beta.threads.runs.create(
        thread_id='thread_Jj24BQoMxDGWvoZw4F7LoeFj',
        assistant_id='asst_hIVzqiJ1sFTkHUssKZrEQymt',  # Ensure 'assistant_response.id' is defined earlier in your script
        instructions="Please address the user as Jane Doe. The user has a premium account."
        )
    print("Run creation response:", run_response)
    
    attempts = 0
    while attempts < 10:  # Limit the number of attempts to check the run's status
        run_status = client.beta.threads.runs.retrieve(
            thread_id=thread_response.id,
            run_id=run_response.id
        )
        print("Run status:", run_status.status)  # Print the current status for debugging
        if run_status.status == "completed":
            print("Run completed.")
            break
        elif run_status.status == "failed":
            print("Run failed.")
            break
        time.sleep(10)  # Poll every 10 seconds
        attempts += 1
        
    if attempts >= 10:
        print("Run did not complete within the expected time.")
    
    messages_response = client.beta.threads.messages.list(
        thread_id=thread_response.id
        )
    print("Messages list response:", messages_response)
    
except Exception as e:
    print("An error occurred:", e)
