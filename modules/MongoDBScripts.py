from pymongo import MongoClient
import pandas as pd

MONGO_CLIENT_ADMIN_URL = ""
MONGO_CLIENT_URL_READONLY = "mongodb+srv://readonlyuser:readonlyuser@cluster0.nrfez.mongodb.net/"

client = MongoClient(MONGO_CLIENT_URL_READONLY)
db = client["FYP"]
fused_features_set_collection = db["fused_features_set"]
cleaned_protein_transcript_liver_set_collection = db[
    "cleaned_protein_transcript_liver_set"
]


def create_appropriate_indexes(columns, collection):
    if "protein_id" in columns and "mrna_id" in columns:
        collection.create_index([("protein_id", 1), ("mrna_id", 1)], unique=True)
    elif "protein_id" in columns:
        collection.create_index("protein_id", unique=True)
    elif "mrna_id" in columns:
        collection.create_index("mrna_id", unique=True)


def csv_to_mongodb(csv_file, collection):
    df = pd.read_csv(csv_file)

    create_appropriate_indexes(df.columns, collection)

    records = df.to_dict("records")

    try:
        collection.insert_many(records, ordered=False)
    except Exception as e:
        print(f"Some documents may be duplicates: {e}")


def fetch_mongodb_to_dataframe(collection, query=None, projection=None):
    if query is None:
        query = {}

    if projection is None:
        cursor = collection.find(query)
    else:
        cursor = collection.find(query, projection)

    data = list(cursor)

    df = pd.DataFrame(data)

    if "_id" in df.columns:
        df = df.drop("_id", axis=1)

    return df
