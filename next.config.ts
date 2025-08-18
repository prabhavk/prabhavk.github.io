// next.config.ts
const isDev = process.env.NODE_ENV !== "production";

const nextConfig = {
  // Do NOT set `output: "export"` on Vercel (you need API routes).
  // output: "export",

  // Your custom key (if you use it elsewhere)
  allowedDevOrigins: ["http://192.168.29.107:3000"],

  // You can add other experimental flags here as needed.
  ...(isDev ? { experimental: {} } : {}),

  async headers() {
    return [
      {
        // Serve WASM & its loader with appropriate caching and isolation
        source: "/wasm/:path*",
        headers: isDev
          ? [
              { key: "Cache-Control", value: "no-store" },
              { key: "Cross-Origin-Opener-Policy", value: "same-origin" },
              { key: "Cross-Origin-Embedder-Policy", value: "require-corp" },
            ]
          : [
              {
                key: "Cache-Control",
                value: "public, max-age=31536000, immutable",
              },
              { key: "Cross-Origin-Opener-Policy", value: "same-origin" },
              { key: "Cross-Origin-Embedder-Policy", value: "require-corp" },
            ],
      },
    ];
  },
};

export default nextConfig;
