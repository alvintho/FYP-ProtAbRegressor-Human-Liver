import pandas as pd
import aiohttp
import asyncio
import time
import math
from collections import defaultdict
from typing import List, Dict

API_BATCH_LIMIT = 50  # Ensembl API limit
CONCURRENT_REQUESTS = 3  # Limit concurrent requests


async def fetch_data(session: aiohttp.ClientSession, url: str, json_data: dict) -> dict:
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    async with session.post(url, headers=headers, json=json_data) as response:
        response.raise_for_status()
        return await response.json()


async def process_batch(
    session: aiohttp.ClientSession, batch: List[str], batch_num: int, total_batches: int
) -> List[Dict]:
    base_url = "https://rest.ensembl.org"
    results = []

    try:
        await asyncio.sleep(1)  # Rate limiting

        lookup_task = fetch_data(session, f"{base_url}/lookup/id", {"ids": batch})
        protein_seq_task = fetch_data(
            session, f"{base_url}/sequence/id", {"ids": batch}
        )

        protein_data, protein_sequence_data = await asyncio.gather(
            lookup_task, protein_seq_task
        )

        transcript_ids = [
            protein_data[pid].get("Parent")
            for pid in batch
            if protein_data.get(pid)
            and any(obj["id"] == pid for obj in protein_sequence_data)
        ]

        mrna_data = []
        mrna_sequence_data = []
        if transcript_ids:
            await asyncio.sleep(1)

            mrna_data = await fetch_data(
                session, f"{base_url}/lookup/id", {"ids": transcript_ids}
            )

            await asyncio.sleep(0.5)

            mrna_sequence_data = await fetch_data(
                session,
                f"{base_url}/sequence/id",
                {"ids": transcript_ids, "type": "cdna"},
            )

            await asyncio.sleep(0.5)

            cds_sequence_data = await fetch_data(
                session,
                f"{base_url}/sequence/id",
                {"ids": transcript_ids, "type": "cds"},
            )

        gene_ids = [
            mrna_data[tid].get("Parent")
            for tid in transcript_ids
            if mrna_data.get(tid)
            and any(obj["id"] == tid for obj in mrna_sequence_data)
        ]

        if gene_ids:
            await asyncio.sleep(1)

            gene_data = await fetch_data(
                session, f"{base_url}/lookup/id", {"ids": gene_ids}
            )

        for pid in batch:
            if protein_data.get(pid) and any(
                obj["id"] == pid for obj in protein_sequence_data
            ):
                protein_info = protein_data[pid]
                transcript_id = protein_info.get("Parent")

                mrna_info = mrna_data[transcript_id]
                gene_id = mrna_info.get("Parent")

                result = {
                    "protein_id": pid,
                    "mrna_id": transcript_id,
                    "gene_id": gene_id,
                    "protein_sequence": next(
                        (
                            seq["seq"]
                            for seq in protein_sequence_data
                            if seq["id"] == pid
                        ),
                        None,
                    ),
                    "protein_length": next(
                        (
                            len(seq["seq"])
                            for seq in protein_sequence_data
                            if seq["id"] == pid
                        ),
                        None,
                    ),
                    "protein_version": next(
                        (
                            seq["version"]
                            for seq in protein_sequence_data
                            if seq["id"] == pid
                        ),
                        None,
                    ),
                    "cdna": (
                        next(
                            (
                                seq["seq"]
                                for seq in mrna_sequence_data
                                if seq["id"] == transcript_id
                            ),
                            None,
                        )
                        if transcript_id and mrna_sequence_data
                        else None
                    ),
                    "cds": (
                        next(
                            (
                                seq["seq"]
                                for seq in cds_sequence_data
                                if seq["id"] == transcript_id
                            ),
                            None,
                        )
                        if transcript_id and cds_sequence_data
                        else None
                    ),
                    "mrna_length": (
                        next(
                            (
                                len(seq["seq"])
                                for seq in mrna_sequence_data
                                if seq["id"] == transcript_id
                            ),
                            None,
                        )
                        if transcript_id and mrna_sequence_data
                        else None
                    ),
                    "mrna_version": (
                        next(
                            (
                                seq["version"]
                                for seq in mrna_sequence_data
                                if seq["id"] == transcript_id
                            ),
                            None,
                        )
                        if transcript_id and mrna_sequence_data
                        else None
                    ),
                    "is_canonical": (
                        mrna_data.get(transcript_id, {}).get("is_canonical")
                        if transcript_id
                        else None
                    ),
                    "canonical_transcript": (
                        gene_data.get(gene_id, {}).get("canonical_transcript")
                        if gene_id
                        else None
                    ),
                    "gene_version": (
                        gene_data.get(gene_id, {}).get("version") if gene_id else None
                    ),
                }
                results.append(result)

        print(
            f"\rProcessing batch {batch_num}/{total_batches} ({len(batch)} proteins)",
            end="",
        )

    except Exception as e:
        print(f"\nError in batch {batch_num}: {str(e)}")

    return results


async def batch_lookup_proteins(protein_ids: List[str]) -> List[Dict]:
    all_results = []
    total_batches = math.ceil(len(protein_ids) / API_BATCH_LIMIT)

    async with aiohttp.ClientSession() as session:
        semaphore = asyncio.Semaphore(CONCURRENT_REQUESTS)
        tasks = []

        for i in range(0, len(protein_ids), API_BATCH_LIMIT):
            batch = protein_ids[i : i + API_BATCH_LIMIT]
            batch_num = (i // API_BATCH_LIMIT) + 1

            async def process_with_semaphore(batch, batch_num):
                async with semaphore:
                    return await process_batch(session, batch, batch_num, total_batches)

            tasks.append(process_with_semaphore(batch, batch_num))

        batch_results = await asyncio.gather(*tasks)

        for results in batch_results:
            all_results.extend(results)

    return all_results


def process_abundance_file(filename: str, max_lines: int = 20000) -> Dict:
    abundance_data = defaultdict(dict)
    line_count = 0

    with open(filename, "r") as f:
        for line in f:
            if line_count >= max_lines:
                break

            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) == 2:
                protein_id = parts[0].split(".")[-1]
                abundance_data[protein_id]["abundance"] = float(parts[1])

            line_count += 1

    return abundance_data


async def main():
    start_time = time.time()

    abundance_data = process_abundance_file(
        "/Users/alvinthosatria/Documents/FYP/FYP/PaxDB/PaxDB_liver_dataset.txt"
    )
    protein_ids = list(abundance_data.keys())

    print(
        f"Processing {len(protein_ids)} proteins in {math.ceil(len(protein_ids)/API_BATCH_LIMIT)} batches"
    )
    results = await batch_lookup_proteins(protein_ids)

    for result in results:
        pid = result["protein_id"]
        if pid in abundance_data:
            result["protein_abundance"] = abundance_data[pid]["abundance"]

    df = pd.DataFrame(results)
    df.to_csv("./liver_paxdb_dataset.csv", index=False)

    print(f"\nTotal execution time: {time.time() - start_time:.2f} seconds")
    print(f"Processed {len(results)} proteins")

    return df


if __name__ == "__main__":
    df = asyncio.run(main())
    print("\nFirst few rows of the results:")
    print(df.head())
