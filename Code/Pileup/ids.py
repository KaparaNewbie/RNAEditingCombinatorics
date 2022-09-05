from math import ceil


# todo instead of skipping `na`, omit `n` s.t.:
# ("a", "z") -> ("a", "m"), ("o", "z"),
# ("A", "Z") -> ("A", "M"), ("O", "Z")
ASCII_RANGES = (
    ("0", "9"),
    ("a", "z"),
    ("A", "Z"),
)


def define_compression_alphabet(ascii_ranges):
    complete_alphabet = []
    for first_letter, last_letter in ascii_ranges:
        local_range_letters = [
            chr(x) for x in range(ord(first_letter), ord(last_letter) + 1)
        ]
        complete_alphabet.extend(local_range_letters)
    return complete_alphabet


def encode(num, alphabet):
    """Encode a positive number into Base X and return the string.

    Source: https://stackoverflow.com/a/1119769/10249633

    Arguments:
    - `num`: The number to encode
    - `alphabet`: The alphabet to use for encoding
    """
    if num == 0:
        return alphabet[0]
    arr = []
    arr_append = arr.append  # Extract bound-method for faster access.
    _divmod = divmod  # Access to locals is faster.
    base = len(alphabet)
    while num:
        num, rem = _divmod(num, base)
        arr_append(alphabet[rem])
    arr.reverse()
    return "".join(arr)


def serialize_compressed_ids(n_ids, ascii_ranges=ASCII_RANGES):
    # compute
    alphabet = define_compression_alphabet(ascii_ranges)
    compressed_ids = []
    x = 0
    while len(compressed_ids) < n_ids:
        compressed_id = encode(x, alphabet)
        # don't create ids such as "na", "NA", "NaN", etc. - they might appear to pandas as missing values
        if "na" not in compressed_id.lower():
            compressed_ids.append(compressed_id)
        x += 1
    # tests
    # (1) uniqueness of ids
    if len(set(compressed_ids)) != len(compressed_ids):
        raise Exception("Encoded ids are not unique")
    # (2) memory efficiency
    max_chars_needed = (
        ceil(len(compressed_ids) ** (1 / len(alphabet))) + 2
    )  # + 1 (instead of + 2) could have been enough if ids with `na` weren't excluded
    max_chars_got = max([len(x) for x in compressed_ids])
    if max_chars_got > max_chars_needed:
        raise Exception(
            f"Encoded ids are not efficient: {max_chars_needed = }, {max_chars_got = }"
        )
    # return list of ids
    return compressed_ids


def compress_ids(old_ids, ascii_ranges=ASCII_RANGES):
    # compute
    n_ids = len(old_ids)
    new_ids = serialize_compressed_ids(n_ids, ascii_ranges)
    # return a mapping from old- to new & compressed ids
    return {old_id: new_id for old_id, new_id in zip(old_ids, new_ids)}
