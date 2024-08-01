FROM continuumio/miniconda3

WORKDIR /code

RUN conda install conda-forge::rdkit

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY . /code

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]

