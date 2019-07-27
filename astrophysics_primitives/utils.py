def sex_to_dec(d, m, s):
    return d + m / 60. + s / 3600.

def time_to_dec(h, m, s):
    return sex_to_dec(h * 15., m, s)

def dec_to_sex(d):
    d = int(d)
    m = (d % 1) * 60
    s = (m % 1) * 60

    return d, int(m), int(s)

def parse_sexagesimal(s, delimiter = ' '):
    s = s.split(delimiter)

# Implement some time methods here
