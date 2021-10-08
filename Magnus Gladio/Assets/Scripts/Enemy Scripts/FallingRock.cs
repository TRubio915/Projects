using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FallingRock : MonoBehaviour
{
    public Rigidbody2D rb; //RB of rock
    public float damage = 25; //Damage rock does to player
    public Transform hitPoint; //hit point
    public float hitRadius; //Hit radius
    public LayerMask whatIsPlayer; //Player layer mask

    [Header("Touch Hit box")]
    public Vector2 knockBackSpeed;
    public Vector2 knockBackSpeedAttack1;
    public Vector2 knockBackSpeedAttack2;
    public Vector2 knockBackSpeedAttack3;
    public Vector2 touchDamageBotLeft;
    public Vector2 touchDamageTopRight;

    public float
        touchDamageWidth,
        touchDamageHeight;

    public Transform touchDamageCheck;

    private float[] attackDetails;
    private bool enableTouchDamage = false;

    // Start is called before the first frame update
    void Start()
    {
        attackDetails = new float[2];
        EnableOverlap();
    }

    private void Update()
    {
        touchDamageBotLeft.Set(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        touchDamageTopRight.Set(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);

        Collider2D hit = Physics2D.OverlapArea(touchDamageBotLeft, touchDamageTopRight, whatIsPlayer);

        if (hit != null && enableTouchDamage)
        {
            hit.GetComponent<Rigidbody2D>().velocity = Vector2.zero;
            attackDetails[0] = damage / 2;
            attackDetails[1] = transform.position.x;
            hit.SendMessage("Damage", attackDetails);

        }
        /*if (GetComponent<Renderer>().isVisible == false)
        {
            Destroy(gameObject);
        }*/
    }

    private void EnableOverlap()
    {
        enableTouchDamage = true;
    }

    public void OnDrawGizmos()
    {
        Gizmos.color = Color.black;
        Vector2 botLeft = new Vector2(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        Vector2 botRight = new Vector2(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y - touchDamageHeight / 2);
        Vector2 topLeft = new Vector2(touchDamageCheck.position.x - touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);
        Vector2 topRight = new Vector2(touchDamageCheck.position.x + touchDamageWidth / 2, touchDamageCheck.position.y + touchDamageHeight / 2);
        Gizmos.DrawLine(botLeft, botRight);
        Gizmos.DrawLine(botRight, topRight);
        Gizmos.DrawLine(topRight, topLeft);
        Gizmos.DrawLine(topLeft, botLeft);
    }

}
