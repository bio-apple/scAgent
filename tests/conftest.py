import os

# Tests must not call external LLMs even if the developer has API keys exported.
os.environ["OPENAI_API_KEY"] = ""
os.environ.pop("ANTHROPIC_API_KEY", None)
