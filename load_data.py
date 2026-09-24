import pandas as pd
import sqlite3


cd = pd.read_csv("data/cell-count.csv")
connection = sqlite3.connect("cells.db")
connection.execute("PRAGMA foreign_keys = ON;")
cur = connection.cursor()
subject_data = cd[
    ["subject", "project", "condition", "age", "sex", "treatment", "response"]
].drop_duplicates(subset = ["subject"])

sample_data = cd[
    ["sample", "sample_type", "subject", "time_from_treatment_start"]
].drop_duplicates(subset = ["sample"])



cell_columns = [
    "b_cell", "cd8_t_cell","cd4_t_cell", "nk_cell", "monocyte"
]

cell_data = pd.melt(cd, id_vars=["sample"], value_vars=cell_columns, 
var_name = "cell_type", value_name = "cell_count")


cur.execute("DROP TABLE IF EXISTS cell_counts")
cur.execute("DROP TABLE IF EXISTS samples")
cur.execute("DROP TABLE IF EXISTS subjects")
cur.execute(""" CREATE TABLE IF NOT EXISTS subjects (
                subject TEXT PRIMARY KEY,
                project TEXT NOT NULL,
                condition TEXT NOT NULL,
                age INTEGER NOT NULL,
                sex TEXT NOT NULL,
                treatment TEXT NOT NULL,
                response TEXT 
            );""")



cur.execute(""" CREATE TABLE IF NOT EXISTS samples (
                sample TEXT PRIMARY KEY,
                sample_type TEXT NOT NULL,
                subject TEXT NOT NULL,
                time_from_treatment_start INTEGER NOT NULL,

                FOREIGN KEY(subject)
                    REFERENCES subjects(subject)

            );""")


cur.execute(""" CREATE TABLE IF NOT EXISTS cell_counts (
                cell_type TEXT NOT NULL,
                sample TEXT NOT NULL,
                cell_count INTEGER NOT NULL,

                PRIMARY KEY (sample, cell_type),

                FOREIGN KEY(sample)
                    REFERENCES samples(sample)
                
            );""")



subject_data.to_sql(
    "subjects",
    connection,
    if_exists="append",
    index=False
)

sample_data.to_sql(
    "samples",
    connection,
    if_exists="append",
    index=False
)

cell_data.to_sql(
    "cell_counts",
    connection,
    if_exists="append",
    index=False
)




connection.commit()
connection.close()



#Part 2

connection = sqlite3.connect("cells.db")
cur = connection.cursor()
connection.execute("PRAGMA foreign_keys = ON;")

p2query = """
    SELECT
        sample,
        SUM(cell_count) OVER(PARTITION BY sample) AS total_count,
        cell_type AS population,
        cell_count AS count,
        ROUND((cell_count * 100.0) / SUM(cell_count) OVER(PARTITION BY sample), 2) AS percentage 
        FROM cell_counts
        ORDER BY sample"""

summary_df = pd.read_sql_query(p2query, connection)

print(summary_df.head(10))

