import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  output: "export",

  experimental: {
    // Explicitly allow LAN dev access (update IP if needed)
    allowedDevOrigins: ["http://192.168.29.107:3000"],
  },
};

export default nextConfig;
