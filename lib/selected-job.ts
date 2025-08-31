// lib/selected-job.ts
export const SELECTED_JOB_KEY = "emtr:selectedJobId";

export function saveSelectedJob(jobId: string) {
  try { localStorage.setItem(SELECTED_JOB_KEY, jobId); } catch {}
}

export function loadSelectedJob(): string {
  try { return localStorage.getItem(SELECTED_JOB_KEY) || ""; } catch { return ""; }
}

export function onSelectedJobChange(cb: (jobId: string) => void) {
  const handler = (e: StorageEvent) => {
    if (e.key === SELECTED_JOB_KEY && typeof e.newValue === "string") cb(e.newValue);
  };
  window.addEventListener("storage", handler);
  return () => window.removeEventListener("storage", handler);
}
