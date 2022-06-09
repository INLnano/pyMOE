
def progress_bar(progress, bar_length=20, bar_character='#'):
    """
    Progress bar.
    Writes a progress bar in place in the output
    
    Args:
        progress: value between 0 and 1
        bar_length: number of characters to consider in the bar
    """
    
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
    block = int(round(bar_length * progress))
#     clear_output(wait = True)
    text = "Progress: [{0}] {1:.1f}%".format( bar_character * block + "-" * (bar_length - block), progress * 100)
    if progress <1:
        end = '\r'
    else:
        end = "\n"
    print(text, end=end) 