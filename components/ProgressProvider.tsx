"use client";

import React, {
  createContext,
  useContext,
  useLayoutEffect, // ⬅️ useLayoutEffect
  useRef,
  useState,
  useCallback,     // ⬅️ useCallback
} from "react";

type Status = "idle" | "running" | "done" | "error";

type ProgressCtx = {
  jobId?: string;
  status: Status;
  logs: string[];
  start: (jobId: string) => void;
  append: (line: string) => void;
  clear: () => void;
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
  const attachScroller = useCallback((el: HTMLElement | null) => {
    scrollerRef.current = el ?? null;
  }, []); // ⬅️ stable function identity

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

  // Auto-scroll to bottom whenever logs/status change
  useLayoutEffect(() => {
    const el = scrollerRef.current;
    if (!el) return;

    const nearBottom = el.scrollHeight - el.scrollTop - el.clientHeight < 40;

    if (nearBottom || status === "running") {
      // wait for layout/paint to account for the new line
      requestAnimationFrame(() => {
        const prev = (el.style as any).scrollBehavior;
        (el.style as any).scrollBehavior = status === "running" ? "smooth" : "auto";
        el.scrollTop = el.scrollHeight;
        (el.style as any).scrollBehavior = prev ?? "";
      });
    }
  }, [logs, status]); // ⬅️ runs after layout thanks to useLayoutEffect

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
