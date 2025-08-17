"use client";

import React, { createContext, useContext, useEffect, useRef, useState } from "react";

type Status = "idle" | "running" | "done" | "error";
type ProgressCtx = {
  jobId?: string;
  status: Status;
  logs: string[];
  start: (jobId: string) => void;
  append: (line: string) => void;
  clear: () => void;
};

const Ctx = createContext<ProgressCtx | null>(null);

export function useProgress() {
  const ctx = useContext(Ctx);
  if (!ctx) throw new Error("useProgress must be used within <ProgressProvider>");
  return ctx;
}

export default function ProgressProvider({ children }: { children: React.ReactNode }) {
  const [jobId, setJobId] = useState<string | undefined>(undefined);
  const [status, setStatus] = useState<Status>("idle");
  const [logs, setLogs] = useState<string[]>([]);
  const pollRef = useRef<number | null>(null);

  function clear() {
    setLogs([]);
    setStatus("idle");
    setJobId(undefined);
    if (pollRef.current) window.clearInterval(pollRef.current);
    pollRef.current = null;
  }

  function append(line: string) {
    setLogs(prev => [...prev, line]);
  }

  // Call this when you kick off EMTR
  function start(id: string) {
    setJobId(id);
    setStatus("running");
    setLogs([]);
    // TODO: Replace with your backend stream (SSE/WebSocket) when ready.
    // --- Demo/poll stub (static hosting fallback) ---
    if (pollRef.current) window.clearInterval(pollRef.current);
    let i = 0;
    pollRef.current = window.setInterval(() => {
      i++;
      append(`[${new Date().toLocaleTimeString()}] EMTR: iteration ${i}, ll=-4363.9…`);
      if (i >= 10) {
        setStatus("done");
        if (pollRef.current) window.clearInterval(pollRef.current);
        pollRef.current = null;        
      }
    }, 700);
  }

  // If you later wire a real backend:
  // useEffect(() => {
  //   if (!jobId) return;
  //   const es = new EventSource(`https://YOUR-BACKEND/jobs/${jobId}/stream`);
  //   es.onmessage = (e) => append(e.data);
  //   es.onerror = () => setStatus("error");
  //   return () => es.close();
  // }, [jobId]);

  return (
    <Ctx.Provider value={{ jobId, status, logs, start, append, clear }}>
      {children}
    </Ctx.Provider>
  );
}
