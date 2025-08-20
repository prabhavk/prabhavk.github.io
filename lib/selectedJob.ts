// lib/selectedJob.ts
export const SELECTED_KEY = "emtr:selectedJobId";

/** Get the currently selected job id from localStorage (client-side). */
export function getSelectedJobId(): string | null {
  if (typeof window === "undefined") return null;
  return window.localStorage.getItem(SELECTED_KEY);
}

/** Subscribe to selection changes across pages/components. Returns an unsubscribe fn. */
export function onSelectedJobChanged(cb: (jobId: string) => void): () => void {
  if (typeof window === "undefined") return () => {};
  const handler = (e: Event) => {
    const ce = e as CustomEvent<{ jobId: string }>;
    if (ce.detail?.jobId) cb(ce.detail.jobId);
  };
  window.addEventListener("emtr:selectedJobChanged", handler as EventListener);
  return () => window.removeEventListener("emtr:selectedJobChanged", handler as EventListener);
}
