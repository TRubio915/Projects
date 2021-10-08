using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraController : MonoBehaviour
{


    public float followSpeed = 2f;
    public float yOffSet = 2.5f;
    public float xOffSet = 0.5f;

    [SerializeField]
    private float
        maxY,
        maxX,
        minY,
        minX;

    public Transform target;
       
    public GameObject border;

    // Start is called before the first frame update

    private float currentY, maxPeekY, minPeekY;

    public bool camPeek = false, enableDeadZone = true;

    public float deadZoneWidth = 1f;

    void Start()
    {
        SetupCamera();
    }

    public void SetupCamera()
    {
        target = GameObject.FindGameObjectWithTag("Player").GetComponent<Transform>();

        updateBorder();
    }

    public void updateBorder()
    {
        minX = border.transform.position.x - (border.GetComponent<BoxCollider2D>().bounds.extents.x);
        maxX = border.transform.position.x + (border.GetComponent<BoxCollider2D>().bounds.extents.x);
        minY = border.transform.position.y - (border.GetComponent<BoxCollider2D>().bounds.extents.y);
        maxY = border.transform.position.y + (border.GetComponent<BoxCollider2D>().bounds.extents.y);
    }

    // Update is called once per frame
    void LateUpdate()
    {
        
        if (target != null && camPeek != true && !GameSettings.isGamePaused)
        {
            Vector3 newPosition = target.position;
           
            newPosition.z = -10;
            newPosition.x = Mathf.Clamp(newPosition.x, minX, maxX);
            newPosition.y = Mathf.Clamp(newPosition.y, minY, maxY);
            

            if (enableDeadZone)
            {
                float deltaX = target.transform.position.x - transform.position.x;
                if (Mathf.Abs(deltaX) > deadZoneWidth)
                {
                    if (deltaX > 0)
                    {
                        newPosition.x = target.transform.position.x - deadZoneWidth;
                    }
                    else
                    {
                        newPosition.x = target.transform.position.x + deadZoneWidth;
                    }
                }
                else
                {
                   
                    newPosition.x = target.transform.position.x - deltaX;
                    
                }
            }
            else
            {
                newPosition.x += xOffSet * target.GetComponent<PlayerController>().GetMovementDirectionOfPlayer();
            }

            
            newPosition.y += yOffSet;

            transform.position = Vector3.Slerp(transform.position, newPosition, followSpeed * Time.deltaTime / Time.timeScale);

            currentY = newPosition.y;
            maxPeekY = currentY + 5;
            minPeekY = currentY - 5;
            
        }   
    }

    public void CamPeek(float direction)
    {
        if (direction == -1 )
        {
            if (currentY > minPeekY)
            {
                currentY = Mathf.MoveTowards(currentY, currentY - 1, 10 * Time.deltaTime);
                transform.position = new Vector3(transform.position.x, currentY, transform.position.z);
            }
        }
        else if (direction == 1)
        {
            if (currentY < maxPeekY)
            {
                currentY = Mathf.MoveTowards(currentY, currentY + 1, 10 * Time.deltaTime);
                transform.position = new Vector3(transform.position.x, currentY, transform.position.z);
            }
        }
        else
            camPeek = false;
    }


}
