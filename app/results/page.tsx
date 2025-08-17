export default function ResultsPage() {
  // Static placeholder for now.
  // Later, hydrate this from your API/DB or read a JSON report.
  const example = [
    { replicate: 1, method: "parsimony", logL: -4370.8162 },
    { replicate: 1, method: "dirichlet", logL: -4363.9295 },
    { replicate: 2, method: "parsimony", logL: -4370.1921 },
    { replicate: 2, method: "dirichlet", logL: -4363.9295 },
  ];

  return (
    <main className="max-w-5xl mx-auto p-6">
      <h1 className="text-2xl font-bold mb-4">Results</h1>
      <p className="text-sm text-gray-600 mb-6">
        This is a placeholder. You can wire this to your backend later.
      </p>

      <div className="overflow-x-auto border rounded-lg">
        <table className="min-w-full text-left text-sm">
          <thead className="bg-gray-100">
            <tr>
              <th className="px-4 py-2">Replicate</th>
              <th className="px-4 py-2">Method</th>
              <th className="px-4 py-2">Log-likelihood</th>
            </tr>
          </thead>
          <tbody>
            {example.map((r, i) => (
              <tr key={i} className="border-t">
                <td className="px-4 py-2">{r.replicate}</td>
                <td className="px-4 py-2">{r.method}</td>
                <td className="px-4 py-2">{r.logL}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Space for future charts */}
      <section className="mt-8">
        <h2 className="text-xl font-semibold mb-2">Visualizations</h2>
        <p className="text-gray-600">Violin &amp; flower plots can go here.</p>
      </section>
    </main>
  );
}
