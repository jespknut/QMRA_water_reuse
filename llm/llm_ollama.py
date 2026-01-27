import requests

OLLAMA_URL = "http://localhost:11434/api/generate"
MODEL = "llama3.1:8b"

SYSTEM_PROMPT = (
    "You are a professional water reuse consultant "
    "writing interpretive QMRA risk summaries."
)

def generate_ai_summary(prompt_text, temperature=0.5):
    payload = {
        "model": MODEL,
        "prompt": f"{SYSTEM_PROMPT}\n\n{prompt_text}",
        "temperature": temperature,
        "stream": False
    }

    r = requests.post(OLLAMA_URL, json=payload, timeout=120)
    r.raise_for_status()

    return r.json()["response"].strip()
