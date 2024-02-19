from datetime import datetime


def final_words():
    final_words_message = "A star shines on the hour of our scripting"
    print(f"\n\n\n{datetime.now()}\t\t{final_words_message}\n", end="")


def ic_prefix():
    return f"ic | {datetime.now()} | "
