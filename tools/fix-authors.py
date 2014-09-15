import sys


def main():
    for line in sys.stdin:
        print fix_name(line.strip())


def fix_name(name):
    if name == 'blueraleigh':
        return 'Mike Grundler'
    if name == 'Dan Rabosky':
        return 'Daniel Rabosky'
    elif name == 'pascaltitle':
        return 'Pascal Title'
    elif name == 'josephwb':
        return 'Joseph W. Brown'
    elif name == 'SimonGreenhill':
        return 'Simon Greenhill'
    else:
        return name


main()
