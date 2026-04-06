from datetime import datetime

from icecream import ic


def final_words():
    final_words_message = "A star shines on the hour of our scripting"
    print(f"\n\n\n{datetime.now()}\t\t{final_words_message}\n", end="")


def ic_prefix():
    return f"ic | {datetime.now()} | "


def configure_ic(*, ansi_color: bool) -> None:
    """
    Configure IceCream (ic) output.
    Default is NO ANSI colors to keep redirected logs readable.
    """
    # Always set prefix/context formatting
    # configure icecream to print the time of the print and the context (file, line, function)
    ic.configureOutput(includeContext=True, prefix=ic_prefix)

    # Ensure enabled, then decide about colors
    try:
        ic.enable()
    except Exception:
        pass

    # Preferred path: native color toggle (supported by newer icecream)
    try:
        ic.configureOutput(color=ansi_color)
        return
    except TypeError:
        # Older icecream: no color= kwarg
        pass

    # Fallback: if we can't reliably turn colors off, avoid emitting ic output at all.
    if not ansi_color:
        try:
            ic.disable()
        except Exception:
            pass
    # else: leave enabled (it will probably colorize, but user asked for ansi_color=True anyway)
