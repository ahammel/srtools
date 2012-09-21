from srtools import sr_sam, sr_postgres
import pytest


class PostgresAlignmentTestSetup(object):
    sam_alignment = sr_sam.SamAlignment("test/test_data/test_insertion.sam")

    pq = "postgres@localhost/test_alignment"

    sr_postgres.postgres_dump(sam_alignment, pq)

    postgres_alignment = sr_postgres.PostgresAlignment(pq)


class TestPostgresAlignmentMethods(PostgresAlignmentTestSetup):
    def test_init(self):
        self.postgres_alignment.rewind()
        self.sam_alignment.rewind()

        for pg_read in self.postgres_alignment:
            sam_read = next(self.sam_alignment)
            assert pg_read == sam_read

        with pytest.raises(StopIteration):
            extra = next(self.sam_alignment)
            print(extra)

    def test_head(self):
        assert self.sam_alignment.head() == self.postgres_alignment.head()


class TestPostgresAlignmentFunctions(PostgresAlignmentTestSetup):
    pass
