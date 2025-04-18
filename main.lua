-- Helper function: Multiply a 3x3 matrix by a ModVector3.
local function mat_mul_vector3(m, v)
     local x = m[1][1] * v.x + m[1][2] * v.y + m[1][3] * v.z
     local y = m[2][1] * v.x + m[2][2] * v.y + m[2][3] * v.z
     local z = m[3][1] * v.x + m[3][2] * v.y + m[3][3] * v.z
     return tm.vector3.Create(x, y, z)
end

-- Helper function: Transpose a 3x3 matrix.
local function mat_transpose(m)
     return {
          { m[1][1], m[2][1], m[3][1] },
          { m[1][2], m[2][2], m[3][2] },
          { m[1][3], m[2][3], m[3][3] }
     }
end

-- A simple Newton-Raphson method for numerically solving
-- f(a) = a * cosh((x - h)/a) + k - y = 0.
local function solve_for_a(x, h, y, k, a_init)
     local a = a_init
     for i = 1, 100 do
          local u = (x - h) / a
          local f = a * ((math.exp(u) + math.exp(-u)) / 2) + k - y
          local df = (math.exp(u) + math.exp(-u)) / 2 - ((x - h) / a) * ((math.exp(u) - math.exp(-u)) / 2)
          if math.abs(df) < 1e-6 then break end
          local a_new = a - f / df
          if math.abs(a_new - a) < 1e-6 then
               a = a_new
               break
          end
          a = a_new
     end
     return a
end

--------------------------------------------------------------------------------
-- Function: generate_catenary_vertical
--
-- Description:
--   Generates a 3D catenary curve such that the sag is vertical (i.e. the
--   cable hangs under gravity).
--
-- Method:
--   - Compute the horizontal displacement by projecting the endpoints onto
--     the XZ-plane.
--   - Define a 2D coordinate system with X along the horizontal direction and
--     Y as the vertical (global up).
--   - Generate the 2D catenary using the cosh function.
--   - Map each 2D point back into 3D by moving horizontally from P1 and adding
--     the vertical offset.
--
-- Parameters:
--   P1  : tm.vector3 -- starting endpoint.
--   P2  : tm.vector3 -- ending endpoint.
--   a   : number     -- initial sag parameter.
--
-- Returns:
--   Array of tm.vector3 representing the 3D catenary points.
--------------------------------------------------------------------------------
local function generate_catenary_vertical(P1, P2, a)
     -- Compute the horizontal vector between P1 and P2 by ignoring Y.
     local horizontalVec = tm.vector3.Create(P2.x - P1.x, 0, P2.z - P1.z)
     local horizontalDistance = horizontalVec:Magnitude()
     if horizontalDistance < 1e-6 then
          horizontalVec = tm.vector3.Create(1, 0, 0)
          horizontalDistance = 1
     else
          horizontalVec = horizontalVec / horizontalDistance
     end

     -- Global up vector: (0,1,0)
     local up = tm.vector3.Up()

     -- Define 2D coordinates for the vertical plane.
     -- Set P1 as (0, P1.y) and P2 as (horizontalDistance, P2.y).
     local x1, y1 = 0, P1.y
     local x2, y2 = horizontalDistance, P2.y

     -- Calculate the midpoint and vertical shift for the catenary equation.
     local h_mid = (x1 + x2) / 2
     local k_shift = (y1 + y2) / 2

     -- Solve for the sag parameter using the first endpoint.
     local a_val = solve_for_a(x1, h_mid, y1, k_shift, a)
     if not a_val or a_val ~= a_val or a_val == math.huge then
          a_val = a
     end

     -- Generate 2D points along the catenary curve.
     local numPoints = 100
     local x_vals, y_vals = {}, {}
     for i = 0, numPoints do
          local t = i / numPoints
          local x_val = x1 * (1 - t) + x2 * t
          local y_val = a_val * ((math.exp((x_val - h_mid) / a_val) + math.exp(-(x_val - h_mid) / a_val)) / 2) + k_shift
          table.insert(x_vals, x_val)
          table.insert(y_vals, y_val)
     end

     -- Map the 2D points back to 3D:
     -- The horizontal displacement is along horizontalVec
     -- and the vertical coordinate is the computed y (with P1.y as the base).
     local catenary3D = {}
     for i = 1, #x_vals do
          local horizontalDisp = horizontalVec * x_vals[i]
          local verticalDisp = up * (y_vals[i] - P1.y)
          catenary3D[i] = P1 + horizontalDisp + verticalDisp
     end

     return catenary3D
end

--------------------------------------------------------------------------------
-- Function: generate_catenary_obj
--
-- (this function thickens the centerline to a ribbon mesh
-- and then exports an OBJ string.)
--------------------------------------------------------------------------------
local function generate_catenary_obj(catenary3D, thickness, uv_scale)
     local vertices   = {} -- "v x y z" lines
     local normals    = {} -- "vn x y z" lines
     local uvs        = {} -- "vt u v" lines
     local faces      = {} -- "f ..." lines

     local numSamples = #catenary3D

     for i = 1, numSamples do
          local tangent
          if i == 1 then
               tangent = catenary3D[2] - catenary3D[1]
          elseif i == numSamples then
               tangent = catenary3D[numSamples] - catenary3D[numSamples - 1]
          else
               tangent = catenary3D[i + 1] - catenary3D[i - 1]
          end
          local tMag = tangent.Magnitude()
          if tMag ~= 0 then
               tangent = tangent / tMag
          end

          local up = tm.vector3.Create(0, 1, 0)
          local side = tangent.Cross(up)
          if side.Magnitude() < 1e-6 then
               up = tm.vector3.Create(0, 0, 1)
               side = tangent.Cross(up)
          end
          local sideMag = side.Magnitude()
          if sideMag ~= 0 then
               side = side / sideMag
          end

          local offset      = side * (thickness / 2)
          local leftVertex  = catenary3D[i] + offset
          local rightVertex = catenary3D[i] - offset

          table.insert(vertices, string.format("v %f %f %f", leftVertex.x, leftVertex.y, leftVertex.z))
          table.insert(vertices, string.format("v %f %f %f", rightVertex.x, rightVertex.y, rightVertex.z))

          local normal = side:Cross(tangent)
          local nMag = normal:Magnitude()
          if nMag ~= 0 then
               normal = normal / nMag
          end
          table.insert(normals, string.format("vn %f %f %f", normal.x, normal.y, normal.z))
          table.insert(normals, string.format("vn %f %f %f", normal.x, normal.y, normal.z))

          local u = ((i - 1) / (numSamples - 1)) * uv_scale
          table.insert(uvs, string.format("vt %f %f", u, 1))
          table.insert(uvs, string.format("vt %f %f", u, 0))
     end

     for i = 1, numSamples - 1 do
          local i1 = (2 * i) - 1
          local i2 = (2 * i)
          local i3 = (2 * i) + 1
          local i4 = (2 * i) + 2

          table.insert(faces, string.format("f %d/%d/%d %d/%d/%d %d/%d/%d",
               i1, i1, i1,
               i3, i3, i3,
               i4, i4, i4))
          table.insert(faces, string.format("f %d/%d/%d %d/%d/%d %d/%d/%d",
               i1, i1, i1,
               i4, i4, i4,
               i2, i2, i2))

          -- Duplicate faces with reversed winding order for the backside:
          table.insert(faces, string.format("f %d/%d/%d %d/%d/%d %d/%d/%d",
               i4, i4, i4,
               i3, i3, i3,
               i1, i1, i1))
          table.insert(faces, string.format("f %d/%d/%d %d/%d/%d %d/%d/%d",
               i2, i2, i2,
               i4, i4, i4,
               i1, i1, i1))
     end

     local obj = table.concat(vertices, "\n") .. "\n" ..
         table.concat(uvs, "\n") .. "\n" ..
         table.concat(normals, "\n") .. "\n" ..
         table.concat(faces, "\n")
     return obj
end

--------------------------------------------------------------------------------
-- Example usage with player input callbacks:

local playerData = {}
local objects = {}

local function updatePos1(callbackData)
     playerData[callbackData.playerId].pos1 = tm.players.GetPlayerTransform(callbackData.playerId).GetPosition()
     playerData[callbackData.playerId].objectPos1.GetTransform().SetPosition(playerData[callbackData.playerId].pos1)
end

local function updatePos2(callbackData)
     playerData[callbackData.playerId].pos2 = tm.players.GetPlayerTransform(callbackData.playerId).GetPosition()
     playerData[callbackData.playerId].objectPos2.GetTransform().SetPosition(playerData[callbackData.playerId].pos2)
end

local function updateA(callbackData)
     if tonumber(callbackData.value) ~= nil then
          if tonumber(callbackData.value) < 1e-2 then
               playerData[callbackData.playerId].a = 1e-2 -- Minimum sag to prevent steepness.
               tm.playerUI.AddSubtleMessageForPlayer(callbackData.playerId, "Cable Tool", "A has been clamped at: 0.01",
                    5)
          elseif tonumber(callbackData.value) > 1e3 then
               playerData[callbackData.playerId].a = 1e3 -- Maximum sag to avoid flatness.
               tm.playerUI.AddSubtleMessageForPlayer(callbackData.playerId, "Cable Tool", "A has been clamped at: 1000",
                    5)
          else
               playerData[callbackData.playerId].a = tonumber(callbackData.value)
          end
          --tm.playerUI.SetUIValue(callbackData.playerId, "txtA", playerData[callbackData.playerId].a)
     end
end
function onPlayerJoined(player)
     playerData[player.playerId] = {
          pos1 = tm.vector3.Create(0, 0, 0),
          pos2 = tm.vector3.Create(0, 0, 0),
          objectPos1 = tm.physics.SpawnObject(tm.vector3.Create(0, 0, 0), "PFB_MovePuzzleBall"),
          objectPos2 = tm.physics.SpawnObject(tm.vector3.Create(0, 0, 0), "PFB_MovePuzzleBall"),
          a = 20,
          thickness = 10,
          uv_scale = 1
     }
     playerData[player.playerId].objectPos1.SetIsStatic(true)
     playerData[player.playerId].objectPos2.SetIsStatic(true)
     tm.playerUI.AddUIButton(player.playerId, "btnPos1", "Update Pos1", updatePos1)
     tm.playerUI.AddUIButton(player.playerId, "btnPos2", "Update Pos2", updatePos2)
     tm.playerUI.AddUIText(player.playerId, "txtA", playerData[player.playerId].a, updateA)
     tm.playerUI.AddUIButton(player.playerId, "btnGenerate", "Generate Catenary", function(callbackData)
          local playerId = callbackData.playerId
          local pos1 = playerData[playerId].pos1
          local pos2 = playerData[playerId].pos2
          local a = playerData[playerId].a
          local thickness = playerData[playerId].thickness
          local uv_scale = playerData[playerId].uv_scale

          local catenary_3d = generate_catenary_vertical(pos1, pos2, a)
          local obj_data = generate_catenary_obj(catenary_3d, thickness, uv_scale)
          tm.os.Log(obj_data)
          tm.os.WriteAllText_Dynamic("catenary" .. #objects .. ".obj", obj_data)
          tm.physics.AddMesh("data_dynamic_willNotBeUploadedToWorkshop/catenary" .. #objects .. ".obj",
               "catenary" .. #objects)
          table.insert(objects, "catenary" .. #objects)
          tm.physics.SpawnCustomObjectConcave(pos1, "catenary" .. #objects - 1)
     end)
end

tm.players.OnPlayerJoined.add(onPlayerJoined)
