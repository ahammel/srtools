import srtools
import sr_postgres

class PostgresAlignmentTestSetup(object):
    sam_alignment = srtools.SamAlignment("test/test_insertion.sam")

    pq = "alex@localhost/test/test_alignment.db?"

    sr_postgres.postgres_dump(sam_alignment, pq)

    postgres_alignment = sr_postgres.PostgresAlignment(pq)


class TestPostgresAlignmentMethods(PostgresAlignmentTestSetup):
    def test_init(self):
        assert self.sam_alignment == self.postgres_alignment


class TestPostgresAlignmentFunctions(PostgresAlignmentTestSetup):
    pass
