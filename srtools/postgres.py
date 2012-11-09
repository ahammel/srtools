"""Tools to convert a SAM file to and from a Postgresql database. Mostly to
demonstrate how sam.Alignment inheritance works

"""
from srtools import Alignment
from srtools.sam import Read
import postgresql


class PostgresAlignment(Alignment):
    """An illumina alignment using data stored as a postgres database. The data
    file is a pg locator for the db.

    """
    def read_generator(self):
        with postgresql.open(self.data_file) as database:
            command = "SELECT * FROM reads ORDER BY id;"
            rows = database.prepare(command)
            for row in rows:
                yield Read("\t".join(str(field) for field in row[1:]))

    def head(self):
        with postgresql.open(self.data_file) as db:
            head_tuple = next(iter(db.prepare("SELECT * FROM head;")))
            return head_tuple[0]


def sql_insert_command(read, table_name, id_number):
    """Returns an SQL command that will insert the given sam.Read into
    the given table. The table is assumed to have the fields specified by the
    postgres_dump function, to which this function is a helper.

    """
    values = [("id", id_number),
              ("qname", read.qname),
              ("flag", read.flag),
              ("rname", read.rname),
              ("pos", read.pos),
              ("mapq", read.mapq),
              ("cigar", str(read.cigar)),
              ("rnext", read.rnext),
              ("pnext", read.pnext),
              ("tlen", read.tlen),
              ("seq", read.seq),
              ("qual", read.qual)]

    field_list = []
    value_list = []
    for field, value in values:
        field_list.append(field)
        if isinstance(value, int):
            value_list.append(str(value))
        else:
            value_list.append("'" + str(value) + "'")

    field_list.append("tags")
    value_list.append("'" + " ".join([str(x) for x in read.tags]) + "'")

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
    with postgresql.open(pq_locator) as database:
        database.execute("DROP TABLE IF EXISTS reads;")
        database.execute("CREATE TABLE reads ( "
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
            database.execute(command)
            id_number += 1

        database.execute("DROP TABLE IF EXISTS head;"
                   "CREATE TABLE head (head  text);")

        head_command = "INSERT INTO head (head) VALUES ('"
        head_command += alignment.head()
        head_command += "');"
        database.execute(head_command)
