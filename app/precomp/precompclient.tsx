// app/precomp/PrecompClient.tsx
"use client";

import { useRouter } from "next/navigation";

type JobRow = {
  job_id: string;
  created_at: string;
  status: string | null;
};

export default function PrecompClient({ jobs }: { jobs: JobRow[] }) {
  const router = useRouter();

  const selectJob = (jobId: string) => {
    try {
      localStorage.setItem("emtr:selectedJobId", jobId);
      // Optional: nudge other tabs; many browsers don’t fire 'storage' in the same tab.
      try {
        window.dispatchEvent(new StorageEvent("storage", { key: "emtr:selectedJobId", newValue: jobId }));
      } catch { /* ignore */ }
    } finally {
      router.push(`/trees?job=${encodeURIComponent(jobId)}`);
    }
  };

  return (
    <div className="mx-auto max-w-5xl p-6 space-y-4">
      <h1 className="text-2xl font-semibold">Precomputed Results</h1>

      {jobs.length === 0 ? (
        <p className="text-gray-300">No jobs found.</p>
      ) : (
        <table className="w-full text-black text-sm border">
          <thead className="bg-gray-200">
            <tr>
              <th className="p-2 text-left">Job ID</th>
              <th className="p-2 text-left">Created</th>
              <th className="p-2 text-left">Status</th>
              <th className="p-2 text-left">Action</th>
            </tr>
          </thead>
          <tbody>
            {jobs.map((j) => (
              <tr key={j.job_id} className="border-t">
                <td className="p-2 font-mono">{j.job_id}</td>
                <td className="p-2">{new Date(j.created_at).toLocaleString()}</td>
                <td className="p-2">{j.status ?? "—"}</td>
                <td className="p-2">
                  <button
                    onClick={() => selectJob(j.job_id)}
                    className="px-3 py-1.5 rounded bg-black text-white hover:bg-gray-900"
                  >
                    Use this job
                  </button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </div>
  );
}
