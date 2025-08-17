// app/page.tsx
"use client";

import { useState } from "react";

export default function Home() {
  const [formData, setFormData] = useState({
    sequenceFile: "",
    topologyFile: "",
    convergenceThreshold: "",
    numRepetitions: "",
    maxIter: "",
  });

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setFormData({
      ...formData,
      [e.target.name]: e.target.value,
    });
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    console.log("Submitted:", formData);
    alert("Form submitted! (Check console for data)");
  };

  return (
    <main className="min-h-screen flex items-center justify-center bg-gray-100 p-6">
      <div className="w-full max-w-lg bg-white shadow-lg rounded-2xl p-8">
        <h1 className="text-2xl font-bold mb-6 text-center text-gray-800">
          EMTR Parameters
        </h1>

        <form onSubmit={handleSubmit} className="space-y-4">
          {/* Sequence File */}
          <div>
            <label className="block text-sm font-medium mb-1">Sequence File</label>
            <input
              type="text"
              name="sequenceFile"
              value={formData.sequenceFile}
              onChange={handleChange}
              className="w-full rounded-lg border-gray-300 shadow-sm focus:border-indigo-500 focus:ring-indigo-500"
              placeholder="Enter sequence file name"
              required
            />
          </div>

          {/* Topology File */}
          <div>
            <label className="block text-sm font-medium mb-1">Topology File</label>
            <input
              type="text"
              name="topologyFile"
              value={formData.topologyFile}
              onChange={handleChange}
              className="w-full rounded-lg border-gray-300 shadow-sm focus:border-indigo-500 focus:ring-indigo-500"
              placeholder="Enter topology file name"
              required
            />
          </div>

          {/* Convergence Threshold */}
          <div>
            <label className="block text-sm font-medium mb-1">Convergence Threshold</label>
            <input
              type="number"
              step="any"
              name="convergenceThreshold"
              value={formData.convergenceThreshold}
              onChange={handleChange}
              className="w-full rounded-lg border-gray-300 shadow-sm focus:border-indigo-500 focus:ring-indigo-500"
              placeholder="e.g. 0.005"
              required
            />
          </div>

          {/* Number of Repetitions */}
          <div>
            <label className="block text-sm font-medium mb-1">Number of Repetitions</label>
            <input
              type="number"
              name="numRepetitions"
              value={formData.numRepetitions}
              onChange={handleChange}
              className="w-full rounded-lg border-gray-300 shadow-sm focus:border-indigo-500 focus:ring-indigo-500"
              placeholder="e.g. 100"
              required
            />
          </div>

          {/* Maximum Iterations */}
          <div>
            <label className="block text-sm font-medium mb-1">Max Iterations</label>
            <input
              type="number"
              name="maxIter"
              value={formData.maxIter}
              onChange={handleChange}
              className="w-full rounded-lg border-gray-300 shadow-sm focus:border-indigo-500 focus:ring-indigo-500"
              placeholder="e.g. 1000"
              required
            />
          </div>

          <button
            type="submit"
            className="w-full bg-indigo-600 text-white py-2 px-4 rounded-lg hover:bg-indigo-700 transition"
          >
            Run EMTR
          </button>
        </form>
      </div>
    </main>
  );
}
