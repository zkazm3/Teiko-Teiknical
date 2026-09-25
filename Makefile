.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py

dashboard:
	python -m streamlit run dashboard/app.py
