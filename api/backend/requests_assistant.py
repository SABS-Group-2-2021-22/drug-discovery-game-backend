from flask import Flask, request, jsonify
import time
from dotenv import load_dotenv
import os
from flask_cors import CORS
from openai import OpenAI


app = Flask(__name__)
CORS(app, resources={r"/api/chat": {"origins": "*"}})
load_dotenv()
client = OpenAI()

@app.route('/api/chat', methods=['POST'])
def chat():
    data = request.json
    prompt = data.get('prompt')  # gets prompt from the user from frontend
    assistant_id = os.getenv("assistant_id")  # Make sure the environment variable name matches
    thread_response = client.beta.threads.create()
    
    message_response = client.beta.threads.messages.create(
        thread_id=thread_response.id,
        role="user",
        content=prompt
    )
    
    run_response = client.beta.threads.runs.create(
        thread_id=thread_response.id,
        assistant_id=assistant_id,
        instructions="The user has a premium account."
    )
    
    attempts = 0
    answer = "No response received."  # Default answer initialization
    while attempts < 15:
        run_status = client.beta.threads.runs.retrieve(
            thread_id=thread_response.id,
            run_id=run_response.id
        )
        if run_status.status == "completed":
            messages_response = client.beta.threads.messages.list(
                thread_id=thread_response.id,
                answer = str(messages_response)
            )
            break
        elif run_status.status == "failed":
            answer = "The run failed to complete successfully."
            break
        elif attempts >= 14:  # Last attempt
            answer = "Run did not complete within the expected time."
        
        time.sleep(10)
        attempts += 1

    return jsonify({'answer': answer})

if __name__ == '__main__':
    app.run(debug=True, port=8000)
