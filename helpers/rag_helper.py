from langchain_openai import OpenAIEmbeddings
from langchain_experimental.text_splitter import SemanticChunker
from langchain import hub
from langchain_core.documents import Document
from typing_extensions import TypedDict, List
from langchain_core.vectorstores import InMemoryVectorStore
import pandas as pd
from langchain_core.prompts import PromptTemplate
from langchain_community.document_loaders import PubMedLoader
from langgraph.constants import START
from langgraph.graph import StateGraph
import os


def create_document_loader(cell_subtypes, experimental_description):
    docs = PubMedLoader(f"{cell_subtypes}" + f"{experimental_description}").load_and_split(
        SemanticChunker(embeddings=OpenAIEmbeddings(model="text-embedding-3-large", api_key=os.environ.get("KEY"))))
    return docs


def initialize_vector_store():
    vector_store = InMemoryVectorStore(
        embedding=OpenAIEmbeddings(model="text-embedding-3-large", api_key=os.environ.get("KEY")))
    return vector_store


def create_template():
    template = '''
    f"Avoid filler statements, mention specific genes and relevant literature, cite at least 7 sources used during RAG retrievals in the format of in-line citations (Author, Year) and at the end of the paper (with proper formatting), making sure that ALL TITLES, AUTHORS, AND DATES ARE FULLY ACCURATE, focus only on immunology pathways, and give a 3  page, dense-paragraph paper in an IMRAD (introduction, methods, results, and discussion) format.
    {context}

    Question: {question}

    Answer: '''

    prompt = PromptTemplate.from_template(template)

    return prompt


def create_graph(vector_store, prompt, llm):
    class State(TypedDict):
        question: str
        context: List[Document]
        answer: str

    def retrieve(state: State):
        retrieved_docs = vector_store.similarity_search(state["question"])
        return {"context": retrieved_docs}

    def generate(state: State):
        docs_content = "\n\n".join(doc.page_content for doc in state["context"])
        messages = prompt.invoke({"question": state["question"], "context": docs_content})
        response = llm.invoke(messages)
        return {"answer": response.content}

    graph_builder = StateGraph(State).add_sequence([retrieve, generate])
    graph_builder.add_edge(START, "retrieve")
    graph = graph_builder.compile()
    return graph


def run_graph(graph, csv_string, cell_subtypes, experimental_description):
    response = graph.invoke({
        "question": f"Given the following CSV of gene groups, whether they are upregulated or downregulated, as well as corresponding genes: \n\n{csv_string}\n\n, as well as that the cell type in question is \n\n{cell_subtypes}\n\n, given that {experimental_description}, propose new mechanisms for why this causes genes to be upregulated and down regulated using other relevant aIMMUNOLOGY pathways found from the data provided."})
    return response['answer']
