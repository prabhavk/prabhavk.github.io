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
  /** Call this with a DOM node that renders the logs to enable auto-scroll */
  attachScroller: (el: HTMLElement | null) => void;
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

  // element to auto-scroll
  const scrollerRef = useRef<HTMLElement | null>(null);
  const attachScroller = (el: HTMLElement | null) => {
    scrollerRef.current = el;
  };

  function clear() {
    setLogs([]);
    setStatus("idle");
    setJobId(undefined);
    if (pollRef.current) window.clearInterval(pollRef.current);
    pollRef.current = null;
  }

  function append(line: string) {
    setLogs((prev) => [...prev, line]);
  }

  // Auto-scroll to bottom whenever logs change
  useEffect(() => {
    const el = scrollerRef.current;
    if (!el) return;

    // If user is near bottom (or we just want to force-scroll), scroll to bottom
    const nearBottom = el.scrollHeight - el.scrollTop - el.clientHeight < 40;
    if (nearBottom || status === "running") {
      el.scrollTop = el.scrollHeight;
    }
  }, [logs, status]);

  // Call this when you kick off EMTR
  function start(id: string) {
    setJobId(id);
    setStatus("running");
    setLogs([]);

    // Demo stub (remove when wiring to worker events)
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

  return (
    <Ctx.Provider value={{ jobId, status, logs, start, append, clear, attachScroller }}>
      {children}
    </Ctx.Provider>
  );
}
