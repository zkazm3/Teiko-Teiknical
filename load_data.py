import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import stats


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


#Part 3



p3query = """
    SELECT
        samples.sample,
        sample_type,
        ROUND((cell_count * 100.0) / SUM(cell_count) OVER(PARTITION BY samples.sample), 2) AS frequency,
        condition,
        response,
        treatment,
        cell_type
        FROM samples
        

        
        JOIN subjects
            ON samples.subject = subjects.subject


        JOIN cell_counts
            ON samples.sample = cell_counts.sample

        WHERE sample_type = 'PBMC'
        AND treatment = 'miraclib'
        AND condition = 'melanoma'

"""

p3_df = pd.read_sql_query(p3query, connection)

cdplot = sns.boxplot(data=p3_df, x = "cell_type", y = "frequency", hue = "response", hue_order = ["no", "yes"])
plt.title("Relative Immune Cell Frequencies of Melanoma Patients Under Miraclib Treatment", fontsize=10)
plt.xlabel("Immune Cell Population")
plt.ylabel("Relative Frequency (%)")

handles, labels = cdplot.get_legend_handles_labels()

cdplot.legend(handles, ["Nonresponder", "Responder"], title = "Treatment Response")

plt.savefig("testBoxPLot.png")

b_cell_response = p3_df[(p3_df["response"] == "yes") & (p3_df["cell_type"] == "b_cell")]["frequency"]
b_cell_nonresponse = p3_df[(p3_df["response"] == "no") & (p3_df["cell_type"] == "b_cell")]["frequency"]

cd8_response = p3_df[(p3_df["response"] == "yes") & (p3_df["cell_type"] == "cd8_t_cell")]["frequency"]
cd8_nonresponse = p3_df[(p3_df["response"] == "no") & (p3_df["cell_type"] == "cd8_t_cell")]["frequency"]

cd4_response = p3_df[(p3_df["response"] == "yes") & (p3_df["cell_type"] == "cd4_t_cell")]["frequency"]
cd4_nonresponse = p3_df[(p3_df["response"] == "no") & (p3_df["cell_type"] == "cd4_t_cell")]["frequency"]

nk_response = p3_df[(p3_df["response"] == "yes") & (p3_df["cell_type"] == "nk_cell")]["frequency"]
nk_nonresponse = p3_df[(p3_df["response"] == "no") & (p3_df["cell_type"] == "nk_cell")]["frequency"]

mono_response = p3_df[(p3_df["response"] == "yes") & (p3_df["cell_type"] == "monocyte")]["frequency"]
mono_nonresponse = p3_df[(p3_df["response"] == "no") & (p3_df["cell_type"] == "monocyte")]["frequency"]

#t_stat, p_value = stats.ttest_ind(b_cell_nonresponse, b_cell_response, equal_var = False)

#Bonferroni correction since the amount of tests run is 5 (very little)
og_p = 0.05
num_tests = 5
bonferroni_p = og_p/num_tests


cell_responses = { 
    "b_cell": (b_cell_response, b_cell_nonresponse),
    "cd8_t_cell": (cd8_response, cd8_nonresponse),
    "cd4_t_cell": (cd4_response, cd4_nonresponse),
    "nk_cell": (nk_response, nk_nonresponse),
    "monocyte": (mono_response, mono_nonresponse)
}

stat_results =  []

for c in cell_responses:
    t_stat,  p_value = stats.ttest_ind(cell_responses[c][0], cell_responses[c][1], equal_var = False)
    stat_results.append({

        "cell_type": c,
        "t_statistic": t_stat,
        "p_value": p_value,
        "significant": p_value < bonferroni_p
    })

    print(f"\nCell Type: {c}")
    print(f"t-statistic: {t_stat}")
    print(f"p-value: {p_value}")

    if p_value < bonferroni_p :
        print("Statistically significant")
    else:
        print("Not statistically significant")
    

stats_df = pd.DataFrame(stat_results)

#Part 4:


p4query = """
    SELECT
        subjects.subject,
        subjects.project,
        condition,
        treatment,
        age,
        sex,
        response,
        sample_type,
        time_from_treatment_start,
        samples.sample
    FROM subjects

    JOIN samples
        ON subjects.subject = samples.subject

    WHERE sample_type = "PBMC"
        AND condition = "melanoma"
        AND time_from_treatment_start = 0
        AND treatment = "miraclib"

    """


p4_df = pd.read_sql_query(p4query, connection)



project_samples =  p4_df["project"].value_counts()
print(project_samples)

response_counts = (
    p4_df.groupby("response")["subject"].nunique()
)

sex_counts = (
    p4_df.groupby("sex")["subject"].nunique()
)
print(sex_counts)

answer = """
        SELECT
            AVG(cell_counts.cell_count) AS average_b_cell_count

        FROM subjects

        JOIN samples
            ON subjects.subject = samples.subject

        JOIN cell_counts
            ON samples.sample = cell_counts.sample
        
        WHERE condition = "melanoma"
        AND sex = "M"
        AND response = "yes"
        AND time_from_treatment_start  = 0
        AND cell_type = "b_cell"

    
    """

answer_df = pd.read_sql_query(answer, connection)

print(answer_df)







