import pandas as pd
import sqlite3


cd = pd.read_csv("data/cell-count.csv")
connection = sqlite3.connect("cells.db")
connection.execute("PRAGMA foreign_keys = ON;")
cur = connection.cursor()
cell_columns = [
    "b_cell", "cd8_t_cell","cd4_t_cell", "nk_cell", "monocyte"
]

cell_test = pd.melt(cd, id_vars=["sample"], value_vars=cell_columns, var_name = "cell_type", value_name = "cell_count")

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
connection.commit()
connection.close()

#Part 2

