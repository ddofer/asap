def apply_mask(array, mask):
    if len(array) == len(mask):
        return [element for element, flag in zip(array, mask) if flag]
    else:
        raise Exception('Cannot apply a mask of a different length')

def bit_to_bool(bit):
    return bit == '1'

def find_first_index_of(string, character_list):
    for i, c in enumerate(string):
        if c in character_list:
            return i

def format_as_csv_value(value):
    if type(value) == bool:
        if value:
            return '1'
        else:
            return '0'
    if type(value) == float:
        return '%.4f' % value
    else:
        return str(value)

def write_csv_line(csv_writer, line):
    csv_writer.writerow(map(format_as_csv_value, line))
