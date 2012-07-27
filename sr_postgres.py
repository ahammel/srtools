import srtools
import postgresql


class PostgresAlignment(srtools.Alignment):
    def read_generator(self):
        with postgresql.open(self.data_file) as db:
            command = "SELECT * FROM reads ORDER BY id;"
            rows = iter(db.prepare(command))
            yield parse_postgres_read(next(rows))


def sql_insert_command(read, table_name, id_number):
    """Returns an SQL command that will insert the given srtools.Read into
    the given table. The table is assumed to have the fields specified by the
    postgres_dump function, to which this function is a helper.

    """
    values = [("id", id_number),
              ("qname", read.qname),
              ("flag", read.flag),
              ("rname", read.rname),
              ("pos", read.pos),
              ("mapq", read.mapq),
              ("cigar", srtools.print_cigar(read.cigar)),
              ("rnext", read.rnext),
              ("pnext", read.pnext),
              ("tlen", read.tlen),
              ("seq", read.seq),
              ("tags", read.tags)]

    field_list = []
    value_list = [] 
    for f, v in values:
        field_list.append(f)
        if isinstance(v, int):
            value_list.append(str(v))
        else:
            value_list.append("'" + str(v) + "'")

    command = "INSERT INTO "
    command += table_name
    command += " ("
    command += ",".join(field_list)
    command += ")"
    command += " VALUES ("
    command += ",".join(value_list)
    command += ");"

    return command


def postgres_dump(alignment, pq_locator):
    """Dumps an alignment of SAM reads into a Postgres database"""
    with postgresql.open(pq_locator) as db:
        db.execute("DROP TABLE reads;")
        db.execute("CREATE TABLE reads ( "
                   "id          int, "
                   "qname       varchar(80), "
                   "flag        int, "
                   "rname       varchar(80), "
                   "pos         int, "
                   "mapq        int, "
                   "cigar       varchar(80), "
                   "rnext       varchar(80), "
                   "pnext       int, "
                   "tlen        int, "
                   "seq         varchar(200), "
                   "qual        varchar(200), "
                   "tags        text"
                   ");")

        id_number = 1
        for read in alignment:
            command = sql_insert_command(read, "reads", id_number)
            db.execute(command)
            id_number += 1
