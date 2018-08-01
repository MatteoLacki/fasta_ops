import csv


def write2csv(filepath, colnames, rows):
    """Dump a sequence to a file.

    Args:
        filepath (str):  The path where to save the iterated sequence of rows.
        colnames (list): List of column names.
        rows (iterable): Iterable of rows to dump.
    """

    with open(filepath, "w") as f:
        writer = csv.writer(f,
                            delimiter=',',
                            quotechar='"',
                            quoting=csv.QUOTE_MINIMAL)
        writer.writerow(colnames)
        for row in rows:
            writer.writerow(row)
