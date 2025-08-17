"use client";

import { useState } from "react";

export default function Home() {
  const [form, setForm] = useState({
    sequenceFile: "",
    topologyFile: "",
    threshold: "",
    repetitions: "",
    maxIter: "",
  });

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setForm({ ...form, [e.target.name]: e.target.value });
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    console.log("Form submitted:", form);
  };

  return (
    <main className="flex min-h-screen flex-col items-center justify-center p-6 bg-white">
      <h1 className="text-2xl font-bold text-black mb-6">Run EMTR</h1>
      <form
        onSubmit={handleSubmit}
        className="w-full max-w-md space-y-4 p-6 border rounded-lg shadow bg-gray-50"
      >
        <div>
          <label className="block mb-1 text-sm font-medium text-black">
            Sequence File
          </label>
          <input
            type="text"
            name="sequenceFile"
            value={form.sequenceFile}
            onChange={handleChange}
            className="w-full border rounded p-2 text-black"
          />
        </div>

        <div>
          <label className="block mb-1 text-sm font-medium text-black">
            Topology File
          </label>
          <input
            type="text"
            name="topologyFile"
            value={form.topologyFile}
            onChange={handleChange}
            className="w-full border rounded p-2 text-black"
          />
        </div>

        <div>
          <label className="block mb-1 text-sm font-medium text-black">
            Convergence Threshold
          </label>
          <input
            type="number"
            step="any"
            name="threshold"
            value={form.threshold}
            onChange={handleChange}
            className="w-full border rounded p-2 text-black"
          />
        </div>

        <div>
          <label className="block mb-1 text-sm font-medium text-black">
            Number of Repetitions
          </label>
          <input
            type="number"
            name="repetitions"
            value={form.repetitions}
            onChange={handleChange}
            className="w-full border rounded p-2 text-black"
          />
        </div>

        <div>
          <label className="block mb-1 text-sm font-medium text-black">
            Maximum Iterations
          </label>
          <input
            type="number"
            name="maxIter"
            value={form.maxIter}
            onChange={handleChange}
            className="w-full border rounded p-2 text-black"
          />
        </div>

        <button
          type="submit"
          className="w-full bg-blue-600 text-white font-semibold py-2 px-4 rounded hover:bg-blue-700"
        >
          Run EMTR
        </button>
      </form>
    </main>
  );
}
