import streamlit as st
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import stats

st.title("Teiko Immune Cell Analysis Dashboard By Zaid Kazmi")

connection = sqlite3.connect("cells.db")

tab2, tab3, tab4 = st.tabs([

    "Part 2: Data Overview",
    "Part 3: Statistical Analysis",
    "Part 4: Data Subset Analysis"

])

with tab2:
    st.header("Treatment Response Frequencies")
    
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
    st.dataframe(summary_df)

    sample_choice = st.selectbox(

        "Choose a sample:",
        summary_df["sample"].unique()
    )

    chosen_sample = summary_df[
        summary_df["sample"] == sample_choice
    ]

    st.dataframe(chosen_sample)


with tab3:
    st.header("Statistical Analysis With Boxplot")

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

    figure, ax = plt.subplots()

    sns.boxplot(
    data = p3_df,
    x = "cell_type",
    y = "frequency",
    hue = "response",
    hue_order = ["no", "yes"],
    ax = ax
    )

    ax.set_xlabel("Immune Cell Population")
    ax.set_ylabel("Relative Frequency (%)")
    ax.set_title("Relative Immune Cell Frequencies of Melanoma Patients Under Miraclib Treatment")

    st.pyplot(figure)

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

    st.subheader("Statistical Test Results")

    stats_df = pd.read_sql_query(
        "SELECT * FROM part3_statistics", connection
    )

    st.dataframe(stats_df)
    st.caption("Two-sample T-tests with Bonferri Correction "
        "(alpha = 0.05 / 5 = 0.01)."
        )

    significant_cells = stats_df[stats_df["significant"] == 1]["cell_type"].tolist()

    st.write("Cell populations with stastically significant differences:",
    ",".join(significant_cells))


with tab4:
    st.header("Baseline Patient Analysis")

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


    response_counts = (
        p4_df.groupby("response")["subject"].nunique()
    )

    sex_counts = (
        p4_df.groupby("sex")["subject"].nunique()
    )

    st.subheader("Project Samples")
    st.dataframe(project_samples)

    st.subheader("Response Counts")
    st.dataframe(response_counts)

    st.subheader("Sex Counts")
    st.dataframe(sex_counts)

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

    average_b = answer_df["average_b_cell_count"].iloc[0]

    st.metric(

        "Average B-Cell Count",
        f"{average_b:,.2f}"
    )

