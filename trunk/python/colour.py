"""
Colour (colour.py)

A module for colourizing text output to the terminal.
"""

# Default colour (reset to this colour)
DEFAULT_CODE = "\033[0;30;47m"

# Dictionary for translating keywords to colour codes
preset_codes = {"default": DEFAULT_CODE, \
                "reset": DEFAULT_CODE, \
                "debug": "\033[0;33m", \
                "warning": "\033[0;33m", \
                "error": "\033[1;31m"}

attributes = {"reset": 0, \
              "bold": 1, \
              "dim": 2, \
              "underline": 4, \
              "blink": 5, \
              "reverse": 7, \
              "hidden": 8}

fg_colours = {"black": 30, \
              "red": 31, \
              "green": 32, \
              "brown": 33, \
              "blue": 34, \
              "purple": 35, \
              "cyan": 36, \
              "white": 37}

bg_colours = {"black": 40, \
              "red": 41, \
              "green": 42, \
              "brown": 43, \
              "blue": 44, \
              "purple": 45, \
              "cyan": 46, \
              "white": 47}

# Current colour
current_code = DEFAULT_CODE


def cset(preset=None, fg='black', bg='white', **attr):
    """Set current colour code.
        If a preset colour code is provided other arguments will
        be ignored.
    """
    global current_code, DEFAULT_CODE
    if preset is not None:
        if preset in preset_codes:
            current_code = preset_codes[preset]
        else:
            print "Unrecognized preset color code!"
    else:
        set_attr = []
        for a in attr.keys():
            if (a in attributes):
                if (attr[a]):
                    set_attr.append(str(attributes[a]))
            else:
                print "Unrecognized attribute: %s" % a
        if len(set_attr) == 0:
            set_attr = ['0']
        if fg in fg_colours:
            fg_val = fg_colours[fg]
        else:
            # assume fb is a code
            fg_val = fg
            
        if bg in bg_colours:
            bg_val = bg_colours[bg]
        else:
            # assume bg is a code
            bg_val = bg
            
        current_code = '\033[%s;%s;%sm' % (";".join(set_attr), \
                            fg_val, bg_val)


def creset():
    """
    Reset current colour code to DEFAULT_CODE.
    """
    cset('default')


def cstring(s, *override, **kwoverride):
    """
    Return the string s with appropriate colour codes
    inserted.

    Keyword is optional. It will override the current
    colour code.
    """
    global current_code, DEFAULT_CODE
    
    # Assume there are overrides
    # Eventually we can check...
    temp_code = current_code
    print override
    print kwoverride
    cset(*override, **kwoverride)
    
    coloured_s = current_code + s + DEFAULT_CODE
    
    current_code = temp_code
    return coloured_s


def cprint(s, *override, **kwoverride):
    """
    Print the string s with appropriate colour codes
    inserted.

    Keyword is optional. It will override the current
    colour code.
    """
    print cstring(s, *override, **kwoverride)


def show_dictionary():
    raise NotImplementedError("colours.show_dictionary needs to be implemented")


def show_colours():
    raise NotImplementedError("colours.show_colours needs to be implemented")


def show_status():
    """Display status of colours module.
        Print global variables.
    """
    # Should we set colour to default?
    print "DEFAULT_CODE:", repr(DEFAULT_CODE)
    print "current_code:", repr(current_code)
    # Should we print a list of all keywords that
    # match the current code?
